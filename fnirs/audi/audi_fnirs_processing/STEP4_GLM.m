%% STEP3.m
% GLM for common vs oddball fNIRS responses
%
% Assumptions & Notes:
% - stim(1) is common (marker=1), stim(2) is oddball (marker=2).
% - If event durations are missing, use stimDurationSec.
% - HbO/HbR labels may be in ml.dataTypeLabel (preferred) or ml.dataType (fallback).
% - Optional block regressor is built from 23-trial blocks if marker=21 is missing.
% - GLM can use AR-IRLS (Barker et al.) or OLS.

% 说明：用于对预处理后的 fNIRS（audi 任务）做单被试 GLM：构造 common vs oddball 事件回归量（可选 block），逐通道拟合 HbO/HbR/HbT，输出回归系数与对比（oddball-common）统计结果和图表。
% Flow
% - 配置区：设置 Homer3 路径、被试/任务列表、stim 索引、HRF 形状、漂移模型、GLM 方法、ROI 通道对、绘图开关等。
% - 发现被试：自动扫描 sub-* 目录；逐被试逐任务读取 *_fnirs_preproc.mat。
% - 构造回归量：从 stim(1)/stim(2) 生成事件序列并与 gamma HRF 卷积；可选 block（marker=21 或 23-trial 推断）。
% - 设计矩阵：[common, oddball, (block), intercept, drift]，漂移可选 none/linear/poly/dct。
% - 选择通道：ROI 源-探测对筛选 + mlActAuto 质量掩码；记录 HbO/HbR/HbT 类型。
% - 时间掩码：tIncAuto 和 tIncAutoCh 去除运动伪迹时间点。
% - GLM 拟合：按通道做 OLS 或 AR-IRLS；计算各回归项 t/p，构造对比 oddball-common。
% - 输出：保存 *_glm.mat、*_glm_table.csv，可选 ROI 汇总表与多张 QC/结果图。


clear; clc;

%% ===================== [1] CONFIG =====================
% Optional: Homer3 path (leave empty if already on MATLAB path)
homer3Dir = '';

% Root derived from this script location
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
rootDir = fileparts(thisDir);

% Subjects / tasks
subjList = {};              % empty = auto-discover sub-*
taskList = {'audi'};        % used for input/output naming
preprocSuffixFmt = '_task-%s_fnirs_preproc.mat';

% Output directories (per-subject, set inside loop)
outDir = '';
figDir = '';

% Stim indices (in preproc MAT)
stimIndexCommon = 1;
stimIndexOddball = 2;

% Event parameters (fallbacks)
stimDurationSec = 0.5;

% Mixed design options
includeBlockRegressor = false; %是否包含block，不放block回归
blockTrials = 23;
blockDurationSec = [];      % [] = auto from median ITI
blockStimName = 'block_onset';

% HRF model (simple gamma)
hrfDurationSec = 30; %HRF 模型的时长窗口（秒），决定卷积后的响应持续多长时间，先用常用值
gammaShape = 6;     %HRF 的 gamma 函数形状与尺度参数，控制上升/峰值与衰减速度。
gammaScale = 1;

% Drift terms
driftModel = 'none';         % 'none' | 'linear' | 'poly' | 'dct'   漂移建模方式；dct 表示用 DCT 基（相当于高通/去低频）
%首选 dct：最稳健、最常用。dctCutoffSec 决定去掉多慢的漂移。常见取值：128s（保守）或 60–120s（更强去漂移）。poly：只在数据很短、或者你明确想用低阶多项式时用（如阶数 2–3）。linear：最弱，只适合很短的记录或漂移极轻微的情况。none：除非你已在前处理阶段做了强高通，否则不建议。
driftPolyOrder = 3;         % used when driftModel='poly'   仅在 driftModel='poly' 时生效，多项式漂移的阶数。
dctCutoffSec = 256;         % 常用128used when driftModel='dct' (high-pass cutoff)  DCT 高通的"截止周期"（秒），值越大，去除的低频越少（更保守）

% GLM method 推荐用ar-irls，目前无法拟合，先用ols
glmMethod = 'ols';      % 'ar-irls' | 'ols'

% AR-IRLS settings (Barker et al.)
arIrlsOrderSec = 0.5         % 用1-2不能收敛 AR model order in seconds (Pmax = round(arIrlsOrderSec/dt)) AR-IRLS 的自回归模型阶数对应的时间尺度（秒），脚本里用Pmax=round(arIrlsOrderSec/dt) 转成点数
                             %设太小：自相关残留;设太大：AR 过拟合

% ROI filtering (source-detector pairs)
useRoiFilter = true;
saveRoiSummary = true;
roiDefs = struct();
roiDefs(1).name = 'PFC';
roiDefs(1).pairs = [ ...
    2 1; 2 2; 3 2; 3 3; 1 3; 2 5; 5 2; 3 6; 6 3; 4 5; ...
    5 5; 5 6; 6 6; 6 4; 5 8; 8 6; 7 8; 8 8; 8 7];
roiDefs(2).name = 'Parietal';
roiDefs(2).pairs = [ ...
    14 18; 18 18; 15 19; 19 19];
roiDefs(3).name = 'Temporal';
roiDefs(3).pairs = [ ...
    9 9; 13 9; 13 13; 9 13; 17 13; 17 17; 13 17; ...
    12 12; 12 16; 16 16; 16 12; 16 20; 20 20; 20 16];

% Plot colors
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue
colorHbT = [0.00 0.60 0.00]; % green

% Plot options
saveHrfFigure = true;

% Debug
processOnlyFirst = false;

%% ===================== [2] ADD HOMER3 TO PATH (OPTIONAL) =====================
if ~isempty(homer3Dir)
    if exist(homer3Dir, 'dir') ~= 7
        error('homer3Dir does not exist: %s', homer3Dir);
    end
    addpath(genpath(homer3Dir));
end

fprintf('[INFO] Assumptions: stim(1)=common, stim(2)=oddball, blockReg=%d\n', includeBlockRegressor);

%% ===================== [3] DISCOVER SUBJECTS =====================
if isempty(subjList)
    subDirs = dir(fullfile(rootDir, 'sub-*'));
    if isempty(subDirs)
        error('No sub-* folders found under: %s', rootDir);
    end
    subjList = {subDirs.name};
end

%% ===================== [4] MAIN LOOP =====================
for si = 1:numel(subjList)
    subName = subjList{si};
    subPath = fullfile(rootDir, subName);

    for ti = 1:numel(taskList)
        taskName = taskList{ti};
        inMat = fullfile(subPath, sprintf([subName preprocSuffixFmt], taskName));

        if exist(inMat, 'file') ~= 2
            fprintf('[SKIP] %s %s: preproc MAT not found.\n', subName, taskName);
            continue;
        end

        fprintf('\n==================== SUBJECT: %s | TASK: %s ====================\n', subName, taskName);
        fprintf('Loading: %s\n', inMat);

        outDir = fullfile(subPath, 'glm');
        figDir = fullfile(subPath, 'glm');
        if ~exist(outDir, 'dir'); mkdir(outDir); end
        if ~exist(figDir, 'dir'); mkdir(figDir); end

        S = load(inMat, 'dc', 'stim', 'mlActAuto', 'tIncAuto', 'tIncAutoCh', 'qc');
        if ~isfield(S, 'dc')
            error('Missing "dc" in %s', inMat);
        end
        if ~isfield(S, 'stim')
            error('Missing "stim" in %s', inMat);
        end

        dc = S.dc;
        stim = S.stim;
        if iscell(stim)
            stim = [stim{:}];
        end

        % --- time & data ---
        t = double(dc.time);
        y = double(dc.dataTimeSeries);
        if size(y,1) ~= numel(t)
            error('Time length (%d) does not match data rows (%d).', numel(t), size(y,1));
        end
        nT = numel(t);

        % --- stim summary ---
        nStim = numel(stim);
        fprintf('[INFO] Stim summary (%d entries):\n', nStim);
        stimNames = cell(nStim,1);
        stimCounts = zeros(nStim,1);
        for k = 1:nStim
            if isprop(stim(k), 'name')
                stimNames{k} = char(stim(k).name);
            elseif isfield(stim(k), 'name')
                stimNames{k} = char(stim(k).name);
            else
                stimNames{k} = sprintf('stim(%d)', k);
            end
            stimData = stim(k).data;
            if isempty(stimData)
                nEv = 0;
            else
                nEv = size(stimData, 1);
            end
            stimCounts(k) = nEv;
            fprintf('  - %s: %d events\n', stimNames{k}, nEv);
        end

        %% ===================== [4b] FIXED STIM INDICES =====================
        stimIndexCommonUse = stimIndexCommon;
        stimIndexOddballUse = stimIndexOddball;

        if nStim < max(stimIndexCommonUse, stimIndexOddballUse)
            error('stim has only %d entries; need stim(%d) and stim(%d).', ...
                nStim, stimIndexCommonUse, stimIndexOddballUse);
        end

        %% ===================== [5] BUILD REGRESSORS =====================
        dt = median(diff(t));
        if ~isfinite(dt) || dt <= 0
            error('Invalid time vector; cannot compute sampling interval.');
        end

        % HRF (gamma)
        tHrf = (0:dt:hrfDurationSec)';
        hrf = (tHrf.^(gammaShape-1)) .* exp(-tHrf/gammaScale) ./ ...
              (gamma(gammaShape) * (gammaScale^gammaShape));
        if max(hrf) > 0
            hrf = hrf ./ max(hrf);
        end

        % Common regressor 构造 common sound 的回归量。先从 stim(1) 取出所有 common
        % 的事件信息，用 uCommon 在时间轴上画出一个"方波事件序列”，然后用 HRF 做卷积：xCommon = conv(uCommon, hrf)，
        % 把神经事件转成预期的血流响应。最后一行 xCommon = xCommon(1:nT) 是把卷积结果裁剪回原始长度。
        stimCommon = stim(stimIndexCommonUse).data;
        if isempty(stimCommon)
            error('stim(%d) is empty.', stimIndexCommonUse);
        end
        uCommon = zeros(nT,1);
        for ei = 1:size(stimCommon,1)
            onset = stimCommon(ei,1);
            if size(stimCommon,2) >= 2 && stimCommon(ei,2) > 0
                dur = stimCommon(ei,2);
            else
                dur = stimDurationSec;
            end
            if size(stimCommon,2) >= 3
                amp = stimCommon(ei,3);
            else
                amp = 1;
            end
            idx = t >= onset & t < (onset + dur);
            uCommon(idx) = uCommon(idx) + amp;
        end
        xCommon = conv(uCommon, hrf);
        xCommon = xCommon(1:nT);

        % Oddball regressor 构造 oddball sound 的回归量
        stimOddball = stim(stimIndexOddballUse).data;
        if isempty(stimOddball)
            error('stim(%d) is empty.', stimIndexOddballUse);
        end
        uOddball = zeros(nT,1);
        for ei = 1:size(stimOddball,1)
            onset = stimOddball(ei,1);
            if size(stimOddball,2) >= 2 && stimOddball(ei,2) > 0
                dur = stimOddball(ei,2);
            else
                dur = stimDurationSec;
            end
            if size(stimOddball,2) >= 3
                amp = stimOddball(ei,3);
            else
                amp = 1;
            end
            idx = t >= onset & t < (onset + dur);
            uOddball(idx) = uOddball(idx) + amp;
        end
        xOddball = conv(uOddball, hrf);
        xOddball = xOddball(1:nT);

        %% ===================== [5b] HRF FIGURE =====================
        if saveHrfFigure
            figH = figure('Visible','off');
            subplot(2,1,1);
            plot(tHrf, hrf, 'k', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('HRF (a.u.)');
            title('HRF (gamma)', 'Interpreter','none');

            subplot(2,1,2);
            plot(t, xCommon, 'Color', colorHbO); hold on;
            plot(t, xOddball, 'Color', colorHbR);
            xlabel('Time (s)');
            ylabel('Regressor (a.u.)');
            title('Convolved regressors: common vs oddball', 'Interpreter','none');
            legend({'Common','Oddball'}, 'Location','best');

            figPathH = fullfile(figDir, sprintf('%s_task-%s_step-glm_hrf.png', subName, taskName));
            saveas(figH, figPathH);
            close(figH);
        end

        % Block duration fallback (from merged events) 估计 block 的持续时间，用于 block 回归项
        dataMerged = [stimCommon; stimOddball];
        [~, sortIdx] = sort(dataMerged(:,1));
        dataMerged = dataMerged(sortIdx,:);
        onsetsAll = dataMerged(:,1);
        iti = diff(onsetsAll);
        iti = iti(isfinite(iti));
        if isempty(iti)
            medianITI = 2.0;
        else
            medianITI = median(iti);
        end
        blockDurationSecFallback = blockDurationSec;
        if isempty(blockDurationSecFallback)
            blockDurationSecFallback = blockTrials * medianITI;
        end

        % Optional block regressor 这段是在构造可选的 block 回归量（用于 mixed GLM）
        xBlock = [];
        if includeBlockRegressor
            blockIdx = find(strcmp(stimNames, blockStimName), 1);
            if ~isempty(blockIdx)
                stimBlock = stim(blockIdx).data;
            else
                % build from merged events (23 trials per block)
                blockStartIdx = 1:blockTrials:size(dataMerged,1);
                blockOnsets = onsetsAll(blockStartIdx);
                stimBlock = [blockOnsets, repmat(blockDurationSecFallback, numel(blockOnsets), 1), ones(numel(blockOnsets),1)];
            end

            uBlock = zeros(nT,1);
            for ei = 1:size(stimBlock,1)
                onset = stimBlock(ei,1);
                if size(stimBlock,2) >= 2 && stimBlock(ei,2) > 0
                    dur = stimBlock(ei,2);
                else
                    dur = blockDurationSecFallback;
                end
                if size(stimBlock,2) >= 3
                    amp = stimBlock(ei,3);
                else
                    amp = 1;
                end
                idx = t >= onset & t < (onset + dur);
                uBlock(idx) = uBlock(idx) + amp;
            end
            xBlock = conv(uBlock, hrf);
            xBlock = xBlock(1:nT);
        end

        %% ===================== [6] DESIGN MATRIX =====================
        %构建 GLM 的设计矩阵 X：先放入两个事件回归量：xCommon 和 xOddball。
        % 如果启用 block 回归量，就再加一列 xBlock。再加一列常数项（intercept），用于拟合整体基线偏移。
        % 如果启用线性漂移，就再加一列时间趋势 linear，用于吸收慢漂移。
        % 最终 X 就是用来回归 HbO/HbR 的自变量矩阵。
        X = [xCommon xOddball];
        regNames = {'common','oddball'};

        if includeBlockRegressor
            X = [X xBlock];
            regNames = [regNames {'block'}];
        end

        % intercept and drift
        X = [X ones(nT,1)];
        regNames = [regNames {'intercept'}];

        driftModelUse = lower(strtrim(driftModel));
        %最后按 driftModel 加漂移项
        switch driftModelUse
            case 'none'
                % no drift terms
            case 'linear'
                tTrend = (t - mean(t));
                X = [X tTrend];
                regNames = [regNames {'linear'}];
            case 'poly'
                if driftPolyOrder < 1
                    warning('driftPolyOrder < 1; no polynomial drift added.');
                else
                    t0 = t - mean(t);
                    tNorm = t0 / max(abs(t0));
                    for po = 1:driftPolyOrder
                        X = [X tNorm.^po];
                        regNames = [regNames {sprintf('poly%d', po)}];
                    end
                end
            case 'dct'
                totalDur = (nT - 1) * dt;
                if ~isfinite(totalDur) || totalDur <= 0 || dctCutoffSec <= 0
                    warning('Invalid total duration; DCT drift skipped.');
                else
                    K = floor((2 * totalDur) / dctCutoffSec);
                    if K < 1
                        fprintf('[INFO] DCT drift: cutoff too long; no DCT terms added.\n');
                    else
                        n = (0:nT-1)';
                        dctReg = zeros(nT, K);
                        for k = 1:K
                            dctReg(:,k) = cos(pi * (n + 0.5) * k / nT);
                            regNames = [regNames {sprintf('dct%02d', k)}];
                        end
                        X = [X dctReg];
                    end
                end
            otherwise
                error('Unknown driftModel: %s', driftModel);
        end

        %% ===================== [7] ROI SELECTION =====================
        ml = dc.measurementList;
        nMeasAll = numel(ml);
        src = [ml.sourceIndex]';
        det = [ml.detectorIndex]';

        if isprop(ml(1), 'dataTypeLabel')
            chrom = string({ml.dataTypeLabel})';
        else
            dataType = [ml.dataType]';
            chrom = strings(nMeasAll,1);
            chrom(dataType==1) = "HbO";
            chrom(dataType==2) = "HbR";
            chrom(dataType==3) = "HbT";
            chrom(chrom=="") = "UNK";
        end

        roiLabel = strings(nMeasAll,1);
        roiMask = true(nMeasAll,1);
        if useRoiFilter
            roiMask = false(nMeasAll,1);
            for ri = 1:numel(roiDefs)
                pairs = roiDefs(ri).pairs;
                pairMask = false(nMeasAll,1);
                for pi = 1:size(pairs,1)
                    pairMask = pairMask | (src==pairs(pi,1) & det==pairs(pi,2));
                end
                roiLabel(pairMask) = roiDefs(ri).name;
                roiMask = roiMask | pairMask;
            end

            if ~any(roiMask)
                error('No channels matched ROI definitions. Check source/det pairs.');
            end
        else
            roiLabel(:) = "ALL";
        end

        %% ===================== [7b] CHANNEL QC MASK (mlActAuto) =====================
        qcMaskAll = true(nMeasAll,1);
        qcNote = 'none';
        mlActVec = [];
        if isfield(S, 'mlActAuto') && ~isempty(S.mlActAuto)
            if iscell(S.mlActAuto)
                mlActVec = S.mlActAuto{1};
            else
                mlActVec = S.mlActAuto;
            end
        end

        if ~isempty(mlActVec)
            mlActVec = mlActVec(:) ~= 0;
            if numel(mlActVec) == nMeasAll
                qcMaskAll = mlActVec;
                qcNote = 'mlActAuto applied to dc';
            else
                warning('mlActAuto size mismatch (got %d, expected %d); QC mask disabled.', ...
                    numel(mlActVec), nMeasAll);
            end
        else
            fprintf('[WARN] mlActAuto not found; using all channels.\n');
        end
        fprintf('[INFO] QC mask: keep %d/%d channels (%s)\n', sum(qcMaskAll), nMeasAll, qcNote);

        glmMask = roiMask & qcMaskAll;
        if ~any(glmMask)
            error('No channels left after ROI+QC filtering.');
        end

        y = y(:, glmMask);
        src = src(glmMask);
        det = det(glmMask);
        chrom = chrom(glmMask);
        roiLabel = roiLabel(glmMask);
        nMeas = numel(src);

        %% ===================== [7c] MOTION TIME MASK (tIncAuto) =====================
        tIncVec = [];
        if isfield(S, 'tIncAuto') && ~isempty(S.tIncAuto)
            if iscell(S.tIncAuto)
                tIncVec = S.tIncAuto{1};
            else
                tIncVec = S.tIncAuto;
            end
        end
        if ~isempty(tIncVec)
            tIncVec = double(tIncVec(:) > 0);
            if numel(tIncVec) ~= nT
                if numel(tIncVec) > nT
                    tIncVec = tIncVec(1:nT);
                    warning('tIncAuto longer than time vector; trimmed to %d.', nT);
                else
                    padN = nT - numel(tIncVec);
                    tIncVec = [tIncVec; ones(padN,1)];
                    warning('tIncAuto shorter than time vector; padded to %d.', nT);
                end
            end
            fprintf('[INFO] tIncAuto: keep %.1f%% timepoints.\n', 100 * mean(tIncVec));
        else
            fprintf('[INFO] tIncAuto not found; using all timepoints.\n');
        end

        tIncCh = [];
        if isfield(S, 'tIncAutoCh') && ~isempty(S.tIncAutoCh)
            if iscell(S.tIncAutoCh)
                tIncCh = S.tIncAutoCh{1};
            else
                tIncCh = S.tIncAutoCh;
            end
        end
        if ~isempty(tIncCh)
            if size(tIncCh,1) ~= nT && size(tIncCh,2) == nT
                tIncCh = tIncCh';
            end
            if size(tIncCh,1) ~= nT
                warning('tIncAutoCh time dimension mismatch (got %d, expected %d); ignoring channel-wise mask.', ...
                    size(tIncCh,1), nT);
                tIncCh = [];
            elseif size(tIncCh,2) ~= nMeasAll
                warning('tIncAutoCh channel dimension mismatch (got %d, expected %d); ignoring channel-wise mask.', ...
                    size(tIncCh,2), nMeasAll);
                tIncCh = [];
            else
                tIncCh = double(tIncCh > 0);
                tIncCh = tIncCh(:, glmMask);
                fprintf('[INFO] tIncAutoCh: keep %.1f%% channel-timepoints (after ROI/QC).\n', ...
                    100 * mean(tIncCh(:)));
            end
        end

        %% ===================== [8] GLM =====================
        glmMethodUse = lower(strtrim(glmMethod));
        switch glmMethodUse
            case 'ar-irls'
                % 用 AR-IRLS (Barker et al.) 回归并计算统计量：
                if exist('ar_glm_final', 'file') ~= 2
                    error('ar_glm_final not found. Add Homer3 to path or set homer3Dir.');
                end
                if exist('robust_ar_fit', 'file') ~= 2
                    error('robust_ar_fit not found. Add Homer3 to path or set homer3Dir.');
                end
                if exist('robustfit', 'file') ~= 2
                    error('robustfit not found. Statistics Toolbox is required for AR-IRLS.');
                end
                if exist('tcdf', 'file') ~= 2
                    error('tcdf not found. Statistics Toolbox is required for AR-IRLS.');
                end
                Pmax = max(1, round(arIrlsOrderSec / dt));
                fprintf('[INFO] GLM method: AR-IRLS (Pmax=%d, %.2fs)\n', Pmax, arIrlsOrderSec);
            case 'ols'
                if exist('tcdf', 'file') ~= 2
                    error('tcdf not found. Statistics Toolbox is required for OLS p-values.');
                end
                fprintf('[INFO] GLM method: OLS\n');
            otherwise
                error('Unknown glmMethod: %s', glmMethod);
        end

        p = size(X,2);
        beta = nan(p, nMeas);
        seReg = nan(p, nMeas);
        tReg = nan(p, nMeas);
        pReg = nan(p, nMeas);
        betaC = nan(1, nMeas);
        seC = nan(1, nMeas);
        tC = nan(1, nMeas);
        pC = nan(1, nMeas);
        dof = nan(1, nMeas);

        % Contrast: oddball - common 做 oddball vs common 的对比
        c = zeros(p,1);
        c(strcmp(regNames, 'oddball')) = 1;
        c(strcmp(regNames, 'common')) = -1;

        %下面是对每个通道做 GLM 拟合，并计算统计量：
        for mi = 1:nMeas
            tMask = true(nT,1);
            if ~isempty(tIncVec)
                tMask = tMask & (tIncVec > 0);
            end
            if ~isempty(tIncCh)
                tMask = tMask & (tIncCh(:,mi) > 0);
            end
            tMask = tMask & all(isfinite(X), 2) & isfinite(y(:,mi));

            nKeep = sum(tMask);
            if nKeep <= p
                warning('Not enough valid timepoints for channel %d; skipping GLM.', mi);
                continue;
            end

            Xi = X(tMask, :);
            yi = y(tMask, mi);

            switch glmMethodUse
                case 'ar-irls'
                    [~, beta_i, tstat_i, pval_i, ~, CovB_i, dfe_i] = ar_glm_final(yi, Xi, Pmax);
                    CovB_i = squeeze(CovB_i);
                    if isempty(CovB_i) || any(~isfinite(CovB_i(:)))
                        warning('CovB invalid for channel %d; skipping stats.', mi);
                        continue;
                    end
                case 'ols'
                    dfe_i = size(Xi,1) - size(Xi,2);
                    if dfe_i <= 0
                        warning('Not enough DOF for channel %d; skipping GLM.', mi);
                        continue;
                    end
                    beta_i = Xi \ yi;
                    resid = yi - (Xi * beta_i);
                    s2 = sum(resid.^2) / dfe_i;
                    XtX = Xi' * Xi;
                    CovB_i = s2 * pinv(XtX);
                    se_i = sqrt(diag(CovB_i));
                    tstat_i = beta_i ./ se_i;
                    pval_i = 2 * tcdf(-abs(tstat_i), dfe_i);
                    if any(~isfinite(se_i)) || any(~isfinite(tstat_i))
                        warning('OLS stats invalid for channel %d; skipping stats.', mi);
                        continue;
                    end
                otherwise
                    error('Unknown glmMethod: %s', glmMethod);
            end

            beta(:, mi) = beta_i(:);
            seReg(:, mi) = sqrt(diag(CovB_i));
            tReg(:, mi) = tstat_i(:);
            pReg(:, mi) = pval_i(:);
            dof(mi) = dfe_i;

            betaC(mi) = c' * beta(:, mi);
            cVar = c' * CovB_i * c;
            if cVar > 0
                seC(mi) = sqrt(cVar);
                tC(mi) = betaC(mi) / seC(mi);
                pC(mi) = 2 * tcdf(-abs(tC(mi)), dof(mi));
            end
        end

        %% ===================== [9] OUTPUT TABLE =====================
        rowCell = {};
        for mi = 1:nMeas
            % per regressor
            for ri = 1:numel(regNames)
                rowCell(end+1,:) = {subName, taskName, mi, src(mi), det(mi), char(chrom(mi)), char(roiLabel(mi)), ...
                    regNames{ri}, beta(ri,mi), seReg(ri,mi), tReg(ri,mi), pReg(ri,mi)}; %#ok<AGROW>
            end
            % contrast
            rowCell(end+1,:) = {subName, taskName, mi, src(mi), det(mi), char(chrom(mi)), char(roiLabel(mi)), ...
                'oddball-common', betaC(mi), seC(mi), tC(mi), pC(mi)}; %#ok<AGROW>
        end

        T = cell2table(rowCell, 'VariableNames', ...
            {'subjId','task','measIndex','src','det','chromophore','roi','term','beta','SE','t','p'});

        outMat = fullfile(outDir, sprintf('%s_task-%s_glm.mat', subName, taskName));
        outCsv = fullfile(outDir, sprintf('%s_task-%s_glm_table.csv', subName, taskName));
        save(outMat, 'beta', 'seReg', 'tReg', 'pReg', 'betaC', 'seC', 'tC', 'pC', ...
            'regNames', 'dof', 'X', 'chrom', 'src', 'det', 'roiLabel', 'roiMask', 'glmMask', 'qcMaskAll', '-v7.3');
        writetable(T, outCsv);

        %% ===================== [10] ROI SUMMARY =====================
        roiSummaryTable = table();
        if saveRoiSummary
            if useRoiFilter
                roiList = {roiDefs.name};
            else
                roiList = {'ALL'};
            end
            chromOrder = ["HbO","HbR","HbT"];
            chromList = chromOrder(ismember(chromOrder, unique(chrom)));
            termList = [regNames {'oddball-common'}];
            roiRows = {};
            for ri = 1:numel(roiList)
                roiName = roiList{ri};
                for ci = 1:numel(chromList)
                    chromName = chromList(ci);
                    chanIdx = (roiLabel == roiName) & (chrom == chromName);
                    for ti2 = 1:numel(termList)
                        termName = termList{ti2};
                        if strcmp(termName, 'oddball-common')
                            vals = betaC(chanIdx);
                        else
                            regIdx = find(strcmp(regNames, termName), 1);
                            vals = beta(regIdx, chanIdx);
                        end
                        nChan = sum(chanIdx);
                        if nChan > 0
                            meanBeta = mean(vals, 'omitnan');
                            sdBeta = std(vals, 'omitnan');
                            seBeta = sdBeta / sqrt(nChan);
                        else
                            meanBeta = NaN;
                            sdBeta = NaN;
                            seBeta = NaN;
                        end
                        roiRows(end+1,:) = {subName, taskName, roiName, char(chromName), termName, nChan, meanBeta, sdBeta, seBeta}; %#ok<AGROW>
                    end
                end
            end
            roiSummaryTable = cell2table(roiRows, 'VariableNames', ...
                {'subjId','task','roi','chromophore','term','nChannels','meanBeta','sdBeta','seBeta'});
            outRoiCsv = fullfile(outDir, sprintf('%s_task-%s_glm_roi_table.csv', subName, taskName));
            writetable(roiSummaryTable, outRoiCsv);
            save(outMat, 'roiSummaryTable', '-append');
        end

        %% ===================== [11] FIGURES =====================
        % QC: retained vs rejected channels if available
        if isfield(S, 'mlActAuto') && ~isempty(S.mlActAuto)
            mlAct = S.mlActAuto;
            if iscell(mlAct)
                if numel(mlAct) == 1
                    mlAct = mlAct{1};
                else
                    try
                        mlAct = cell2mat(mlAct(:));
                    catch
                        mlAct = mlAct{1};
                        warning('mlActAuto is a cell array; using first cell for QC plot.');
                    end
                end
            end
            if isnumeric(mlAct) || islogical(mlAct)
                mlActVec = mlAct(:);
                if numel(mlActVec) == nMeasAll
                    mlActVec = mlActVec(roiMask);
                    nActive = sum(mlActVec > 0);
                    nReject = numel(mlActVec) - nActive;
                else
                    warning('mlActAuto size does not match measurement count; skip QC bar plot.');
                    nActive = [];
                    nReject = [];
                end
                fig1 = figure('Visible','off');
                if ~isempty(nActive)
                    bar([nActive nReject]);
                else
                    bar([0 0]);
                end
                set(gca, 'XTickLabel', {'Active','Rejected'});
                ylabel('Measurements');
                title(sprintf('%s %s QC channels', subName, taskName), 'Interpreter','none');
                figPath1 = fullfile(figDir, sprintf('%s_task-%s_step-glm_qc_channels.png', subName, taskName));
                saveas(fig1, figPath1);
                close(fig1);
            else
                warning('mlActAuto is not numeric/logical; skip QC bar plot.');
            end
        end

        % Signal plot: mean HbO/HbR/HbT across channels
        fig2 = figure('Visible','off');
        hboIdx = chrom == "HbO";
        hbrIdx = chrom == "HbR";
        hbtIdx = chrom == "HbT";
        hLine = [];
        legCell = {};
        if any(hboIdx)
            hLine(end+1) = plot(t, mean(y(:,hboIdx), 2, 'omitnan'), 'Color', colorHbO); %#ok<AGROW>
            legCell{end+1} = 'HbO'; %#ok<AGROW>
            hold on;
        end
        if any(hbrIdx)
            hLine(end+1) = plot(t, mean(y(:,hbrIdx), 2, 'omitnan'), 'Color', colorHbR); %#ok<AGROW>
            legCell{end+1} = 'HbR'; %#ok<AGROW>
            hold on;
        end
        if any(hbtIdx)
            hLine(end+1) = plot(t, mean(y(:,hbtIdx), 2, 'omitnan'), 'Color', colorHbT); %#ok<AGROW>
            legCell{end+1} = 'HbT'; %#ok<AGROW>
        end
        xlabel('Time (s)');
        ylabel('Hb (a.u.)');
        title(sprintf('%s %s mean HbO/HbR/HbT', subName, taskName), 'Interpreter','none');
        if ~isempty(hLine)
            legend(hLine, legCell, 'Location','best');
        end
        figPath2 = fullfile(figDir, sprintf('%s_task-%s_step-glm_signal.png', subName, taskName));
        saveas(fig2, figPath2);
        close(fig2);

        % Result plot: contrast distribution
        fig3 = figure('Visible','off');
        if any(hboIdx) || any(hbrIdx) || any(hbtIdx)
            dataPlot = [];
            grp = {};
            if any(hboIdx)
                dataPlot = [dataPlot; betaC(hboIdx)'];
                grp = [grp; repmat({'HbO'}, sum(hboIdx), 1)];
            end
            if any(hbrIdx)
                dataPlot = [dataPlot; betaC(hbrIdx)'];
                grp = [grp; repmat({'HbR'}, sum(hbrIdx), 1)];
            end
            if any(hbtIdx)
                dataPlot = [dataPlot; betaC(hbtIdx)'];
                grp = [grp; repmat({'HbT'}, sum(hbtIdx), 1)];
            end
            boxplot(dataPlot, grp);
            ylabel('Beta (oddball-common)');
        else
            plot(betaC, 'k.');
            ylabel('Beta (oddball-common)');
        end
        title(sprintf('%s %s contrast distribution', subName, taskName), 'Interpreter','none');
        figPath3 = fullfile(figDir, sprintf('%s_task-%s_step-glm_result.png', subName, taskName));
        saveas(fig3, figPath3);
        close(fig3);

        % ROI plot: oddball-common per ROI
        if saveRoiSummary && useRoiFilter
            roiList = {roiDefs.name};
            chromOrder = ["HbO","HbR","HbT"];
            chromList = chromOrder(ismember(chromOrder, unique(chrom)));
            roiMeans = nan(numel(roiList), numel(chromList));
            for ri = 1:numel(roiList)
                for ci = 1:numel(chromList)
                    chanIdx = (roiLabel == roiList{ri}) & (chrom == chromList(ci));
                    if any(chanIdx)
                        roiMeans(ri,ci) = mean(betaC(chanIdx), 'omitnan');
                    end
                end
            end
            fig4 = figure('Visible','off');
            b = bar(categorical(roiList), roiMeans);
            if numel(b) >= 1
                b(1).FaceColor = colorHbO;
            end
            if numel(b) >= 2
                b(2).FaceColor = colorHbR;
            end
            if numel(b) >= 3
                b(3).FaceColor = colorHbT;
            end
            ylabel('Beta (oddball-common)');
            legend(cellstr(chromList), 'Location','best');
            title(sprintf('%s %s ROI contrast', subName, taskName), 'Interpreter','none');
            figPath4 = fullfile(figDir, sprintf('%s_task-%s_step-glm_roi_result.png', subName, taskName));
            saveas(fig4, figPath4);
            close(fig4);
        end

        fprintf('[DONE] %s %s | out=%s\n', subName, taskName, outMat);

        if processOnlyFirst
            break;
        end
    end

    if processOnlyFirst
        break;
    end
end
