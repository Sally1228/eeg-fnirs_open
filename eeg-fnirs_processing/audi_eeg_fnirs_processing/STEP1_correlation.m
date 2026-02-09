% EEG-fNIRS feature correlation (group-level) for the auditory oddball task.
% Assumptions & Notes:
% - EEG features (P300 amplitude/latency + band power) are computed here from epoched .set files.
% - EEG bins: common=1, oddball=2 (update if your ERPLAB bin list differs).
% - fNIRS features are estimated here with a simple GLM on preprocessed dc/stim using a gamma HRF.
% - Analyses include: All stimuli, Common only, Oddball-Common difference.
% - Alignment summary uses STEP2 alignment CSVs (oddball trials only).
% - Short ISI implies strong HRF overlap; interpret "event-level" coupling as GLM-derived.
% - Scatter plots are generated from the computed features.

% 说明：
% - 给出三种条件：All（所有刺激）、Common、OddballCommon（oddball−common）。
% - 组水平相关 + 保存散点图。
% 1) 配置区
% - 定义 EEG 的 .set 文件模式、bin 编号、通道、P300窗口/基线/极性、频段与功率窗口。
% - 定义 fNIRS 的 stim 索引、HRF 形状、漂移模型（DCT/多项式等）、ROI 通道对与 HbO/HbR。
% - 定义相关分析与散点图输出设置。
% 2) 被试发现与路径准备
% - 自动从 eeg/audi/sub-* 和 fnirs/audi/sub-* 取交集作为被试列表。
% - 输出目录：eeg_fnirs_processing/output/STEP1_correlation/。
% 3) EEG特征计算（关键）
% - 对每个被试读取分段 .set。
% - 剔除 EEG.reject 标记的坏试次。
% - 根据 bin 得到三类试次：All / Common / Oddball。
% 在指定通道上计算：
% - P300 幅度与潜伏期：All、Common、Oddball−Common 三个版本。
% - Delta/Theta 功率：同样分 All / Common / Oddball−Common。
% - 结果汇总到 eegFeatures 表。
% 4) fNIRS GLM 特征计算（关键）
% - 直接读取 *_fnirs_preproc.mat，用 dc 和 stim 构造 HRF 回归量：
% - xAll = xCommon + xOdd（所有刺激）
% - xCommon、xOdd（分别）
% - 在 ROI × HbO/HbR 上做 OLS 拟合：
% - 得到 betaAll、betaCommon、betaOdd，再算 betaOdd - betaCommon。
% - 结果汇总到 fnirsFeatures 表，并标注 Analysis=All/Common/OddballCommon。
% 5) 对齐统计（可选）
% - 读取 ..._eeg_fnirs_alignment.csv。
% - 仅用于输出 oddball 的有效/匹配试次数量，做质量检查。
% 6) 相关分析
% - 对每个条件（All/Common/OddballCommon），每个 EEG 特征 × ROI × Hb 类型做 Pearson 相关。
% - 结果输出到 STEP1_group_correlation.csv。
% 7) 散点图
% -基于上述相关结果绘制散点 + 拟合线，保存到 eeg_fnirs_processing/output/STEP1_correlation/figures/。


%% Config
close all; clearvars; clc;

taskName = 'audi';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset patterns (epoched)
eegSetPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set', ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp.set', ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch.set' ...
    };

% EEG bin mapping
eegCommonBin = 1;
eegOddballBin = 2;

% EEG channel selection (ERP + power)
eegChanLabels = {'Pz'};
eegChanFallbackIdx = 53;
eegCombineMode = 'mean'; % 'mean' | 'median'

% P300 settings
p300WindowMs = [300 500];
baselineWindowMs = [-200 0];
p300Measure = 'peak'; % 'peak' | 'mean'
p300Polarity = 'positive'; % 'positive' | 'negative' | 'abs'
baselineMode = 'subtract'; % 'subtract' | 'none'

% EEG power settings
eegPowerWindowMs = [200 500];
eegPowerBaselineWindowMs = [-200 0];
eegPowerBaselineMode = 'none'; % 'none' | 'subtract' | 'db'
eegPowerBandNames = {'Delta', 'Theta'};
eegPowerBandDefs = [1 3; 4 7];

% fNIRS stim settings
stimIndexCommon = 1;
stimIndexOddball = 2;
stimDurationSec = 0.5;
applyFnirsTIncAuto = true;

% fNIRS HRF/GLM settings
hrfDurationSec = 30;
gammaShape = 6;
gammaScale = 1;
driftModel = 'dct'; % 'none' | 'linear' | 'poly' | 'dct'
driftPolyOrder = 3;
dctCutoffSec = 128;

% fNIRS ROI settings
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
chromList = {'HbO', 'HbR'};
minChannels = 1;

% Alignment summary (uses STEP2 output)
useAlignmentSummary = true;
alignmentDirName = 'STEP2_p300_informed_fnirs';

% Correlation settings
analysisSuffix = {'All', 'Common', 'OddballCommon'};
analysisNames = {'All', 'Common', 'OddballCommon'};
minSubjectsForCorr = 3;

% Scatter plot settings
saveScatterPlots = true;
scatterMaxPlots = 60;
scatterFigureFormat = 'png';
scatterDpi = 150;
scatterVisible = 'off'; % 'on' | 'off'
scatterFitLine = true;

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', 'STEP1_correlation');
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end
alignmentRoot = fullfile(scriptDir, 'output', alignmentDirName);

fprintf('[%s] Assumptions: EEG features computed from epoched datasets; fNIRS GLM uses gamma HRF (%s).\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), driftModel);

%% Discover subjects
if isempty(subjList)
    eegSubs = dir(fullfile(eegRoot, 'sub-*'));
    fnirsSubs = dir(fullfile(fnirsRoot, 'sub-*'));
    eegNames = {eegSubs.name};
    fnirsNames = {fnirsSubs.name};
    subjList = intersect(eegNames, fnirsNames, 'stable');
end
if isempty(subjList)
    error('No subjects found under %s and %s', eegRoot, fnirsRoot);
end

%% Compute EEG features (ERP + band power)
if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %#ok<ASGLU>

% Build EEG feature names
featureNames = {};
baseErpNames = {'P300_MeanAmp', 'P300_PeakLatencyMs'};
for bi = 1:numel(baseErpNames)
    for ai = 1:numel(analysisSuffix)
        featureNames{end+1} = sprintf('%s_%s', baseErpNames{bi}, analysisSuffix{ai}); %#ok<AGROW>
    end
end
basePowNames = strcat(eegPowerBandNames, 'Power');
for bi = 1:numel(basePowNames)
    for ai = 1:numel(analysisSuffix)
        featureNames{end+1} = sprintf('%s_%s', basePowNames{bi}, analysisSuffix{ai}); %#ok<AGROW>
    end
end

eegVarNames = [{'SubID'}, featureNames];
eegRows = cell(0, numel(eegVarNames));

for si = 1:numel(subjList)
    subName = subjList{si};
    fprintf('[%s] EEG features: %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, si, numel(subjList));

    eegSubDir = fullfile(eegRoot, subName);
    eegSetPath = '';
    for pi = 1:numel(eegSetPatterns)
        candidate = sprintf(eegSetPatterns{pi}, subName, taskName);
        candidatePath = fullfile(eegSubDir, candidate);
        if exist(candidatePath, 'file') == 2
            eegSetPath = candidatePath;
            break;
        end
    end
    if isempty(eegSetPath)
        fprintf('[%s] WARNING: EEG .set not found for %s. Skipping EEG features.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    [setDir, setName, setExt] = fileparts(eegSetPath);
    EEG = pop_loadset('filename', [setName setExt], 'filepath', setDir);
    if EEG.trials <= 1
        fprintf('[%s] WARNING: EEG dataset is not epoched for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % Remove epochs marked for artifact rejection
    reject = false(1, EEG.trials);
    if isfield(EEG, 'reject') && ~isempty(EEG.reject)
        reject_fields = fieldnames(EEG.reject);
        for f = 1:numel(reject_fields)
            val = EEG.reject.(reject_fields{f});
            if isempty(val) || ~(islogical(val) || isnumeric(val))
                continue;
            end
            if isvector(val) && numel(val) == EEG.trials
                reject = reject | logical(val(:))';
            elseif ismatrix(val) && size(val, 2) == EEG.trials
                reject = reject | any(val, 1);
            elseif ismatrix(val) && size(val, 1) == EEG.trials
                reject = reject | any(val, 2)';
            end
        end
    end
    if any(reject)
        EEG = pop_rejepoch(EEG, reject, 0);
    end

    % Build epoch masks for bins
    commonMask = false(1, EEG.trials);
    oddMask = false(1, EEG.trials);
    if isfield(EEG, 'epoch') && ~isempty(EEG.epoch)
        for e = 1:EEG.trials
            if isfield(EEG.epoch(e), 'eventbini')
                bini = EEG.epoch(e).eventbini;
                if iscell(bini)
                    bini = [bini{:}];
                end
                if any(bini == eegCommonBin)
                    commonMask(e) = true;
                end
                if any(bini == eegOddballBin)
                    oddMask(e) = true;
                end
            end
        end
    end
    allMask = commonMask | oddMask;
    if ~any(allMask)
        fprintf('[%s] WARNING: No common/oddball epochs for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % Channel selection
    chanIdx = [];
    if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
        labels = {EEG.chanlocs.labels};
        for li = 1:numel(eegChanLabels)
            idx = find(strcmpi(labels, eegChanLabels{li}), 1);
            if ~isempty(idx)
                chanIdx = idx;
                break;
            end
        end
    end
    if isempty(chanIdx)
        if eegChanFallbackIdx <= EEG.nbchan
            chanIdx = eegChanFallbackIdx;
        else
            fprintf('[%s] WARNING: Channel index invalid for %s. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
            continue;
        end
    end
    EEG = pop_select(EEG, 'channel', chanIdx);

    data = EEG.data;
    if strcmpi(eegCombineMode, 'median')
        dataChan = squeeze(median(data, 1));
    else
        dataChan = squeeze(mean(data, 1));
    end
    if size(dataChan, 1) ~= numel(EEG.times)
        dataChan = dataChan';
    end
    times = double(EEG.times(:));

    % Baseline correction for ERP
    baseIdx = times >= baselineWindowMs(1) & times <= baselineWindowMs(2);
    if strcmpi(baselineMode, 'subtract') && any(baseIdx)
        baseVal = mean(dataChan(baseIdx, :), 1, 'omitnan');
        dataChan = dataChan - baseVal;
    end

    % P300 time window
    p300Idx = times >= p300WindowMs(1) & times <= p300WindowMs(2);
    if ~any(p300Idx)
        fprintf('[%s] WARNING: P300 window outside data range for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end

    p300AmpAll = NaN; p300LatAll = NaN;
    p300AmpCommon = NaN; p300LatCommon = NaN;
    p300AmpDiff = NaN; p300LatDiff = NaN;

    if any(allMask)
        erpAll = mean(dataChan(:, allMask), 2, 'omitnan');
        seg = erpAll(p300Idx);
        segT = times(p300Idx);
        if ~isempty(seg)
            if strcmpi(p300Measure, 'mean')
                p300AmpAll = mean(seg, 'omitnan');
                p300LatAll = NaN;
            else
                switch lower(p300Polarity)
                    case 'positive'
                        [p300AmpAll, iMax] = max(seg);
                    case 'negative'
                        [p300AmpAll, iMax] = min(seg);
                    otherwise
                        [~, iMax] = max(abs(seg));
                        p300AmpAll = abs(seg(iMax));
                end
                p300LatAll = segT(iMax);
            end
        end
    end

    if any(commonMask)
        erpCommon = mean(dataChan(:, commonMask), 2, 'omitnan');
        seg = erpCommon(p300Idx);
        segT = times(p300Idx);
        if ~isempty(seg)
            if strcmpi(p300Measure, 'mean')
                p300AmpCommon = mean(seg, 'omitnan');
                p300LatCommon = NaN;
            else
                switch lower(p300Polarity)
                    case 'positive'
                        [p300AmpCommon, iMax] = max(seg);
                    case 'negative'
                        [p300AmpCommon, iMax] = min(seg);
                    otherwise
                        [~, iMax] = max(abs(seg));
                        p300AmpCommon = abs(seg(iMax));
                end
                p300LatCommon = segT(iMax);
            end
        end
    else
        erpCommon = [];
    end

    if any(oddMask)
        erpOdd = mean(dataChan(:, oddMask), 2, 'omitnan');
    else
        erpOdd = [];
    end

    if ~isempty(erpOdd) && ~isempty(erpCommon)
        erpDiff = erpOdd - erpCommon;
        seg = erpDiff(p300Idx);
        segT = times(p300Idx);
        if ~isempty(seg)
            if strcmpi(p300Measure, 'mean')
                p300AmpDiff = mean(seg, 'omitnan');
                p300LatDiff = NaN;
            else
                switch lower(p300Polarity)
                    case 'positive'
                        [p300AmpDiff, iMax] = max(seg);
                    case 'negative'
                        [p300AmpDiff, iMax] = min(seg);
                    otherwise
                        [~, iMax] = max(abs(seg));
                        p300AmpDiff = abs(seg(iMax));
                end
                p300LatDiff = segT(iMax);
            end
        end
    end

    % Band power (trial-level, then average by condition)
    timeIdxPower = times >= eegPowerWindowMs(1) & times <= eegPowerWindowMs(2);
    baseIdxPower = times >= eegPowerBaselineWindowMs(1) & times <= eegPowerBaselineWindowMs(2);
    if ~any(timeIdxPower)
        fprintf('[%s] WARNING: Power window outside data range for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end

    numBands = size(eegPowerBandDefs, 1);
    powerAll = nan(1, numBands);
    powerCommon = nan(1, numBands);
    powerDiff = nan(1, numBands);

    for b = 1:numBands
        band = eegPowerBandDefs(b, :);
        if exist('pop_eegfiltnew', 'file') == 2
            EEG_filt = pop_eegfiltnew(EEG, 'locutoff', band(1), 'hicutoff', band(2), 'plotfreqz', 0);
        else
            if exist('eegfilt', 'file') ~= 2
                error('Neither pop_eegfiltnew nor eegfilt is available on the MATLAB path.');
            end
            EEG_filt = EEG;
            data_filt = zeros(size(EEG.data), 'like', EEG.data);
            for tIdx = 1:EEG.trials
                data_filt(:, :, tIdx) = eegfilt(EEG.data(:, :, tIdx), EEG.srate, band(1), band(2));
            end
            EEG_filt.data = data_filt;
        end

        dataPow = EEG_filt.data;
        if strcmpi(eegCombineMode, 'median')
            dataPow = squeeze(median(dataPow, 1));
        else
            dataPow = squeeze(mean(dataPow, 1));
        end
        if size(dataPow, 1) ~= numel(times)
            dataPow = dataPow';
        end

        powTw = mean(dataPow(timeIdxPower, :).^2, 1, 'omitnan');
        powTw = powTw(:)';

        if ~strcmpi(eegPowerBaselineMode, 'none') && any(baseIdxPower)
            powBase = mean(dataPow(baseIdxPower, :).^2, 1, 'omitnan');
            powBase = powBase(:)';
            switch lower(eegPowerBaselineMode)
                case 'subtract'
                    powTw = powTw - powBase;
                case 'db'
                    powTw = 10 * log10((powTw + eps) ./ (powBase + eps));
            end
        end

        if any(allMask)
            powerAll(b) = mean(powTw(allMask), 'omitnan');
        end
        if any(commonMask)
            powerCommon(b) = mean(powTw(commonMask), 'omitnan');
        end
        if any(oddMask) && any(commonMask)
            powerDiff(b) = mean(powTw(oddMask), 'omitnan') - mean(powTw(commonMask), 'omitnan');
        end
    end

    vals = struct();
    vals.P300_MeanAmp_All = p300AmpAll;
    vals.P300_PeakLatencyMs_All = p300LatAll;
    vals.P300_MeanAmp_Common = p300AmpCommon;
    vals.P300_PeakLatencyMs_Common = p300LatCommon;
    vals.P300_MeanAmp_OddballCommon = p300AmpDiff;
    vals.P300_PeakLatencyMs_OddballCommon = p300LatDiff;

    for b = 1:numBands
        baseName = sprintf('%sPower', eegPowerBandNames{b});
        vals.(sprintf('%s_All', baseName)) = powerAll(b);
        vals.(sprintf('%s_Common', baseName)) = powerCommon(b);
        vals.(sprintf('%s_OddballCommon', baseName)) = powerDiff(b);
    end

    row = cell(1, numel(eegVarNames));
    row{1} = subName;
    for vi = 1:numel(featureNames)
        key = featureNames{vi};
        if isfield(vals, key)
            row{vi+1} = vals.(key);
        else
            row{vi+1} = NaN;
        end
    end
    eegRows(end+1, :) = row; %#ok<SAGROW>
end

if isempty(eegRows)
    error('No EEG features computed. Check EEG paths and settings.');
end

eegFeatures = cell2table(eegRows, 'VariableNames', eegVarNames);

%% Compute fNIRS GLM features (ROI-level)
fnirsRows = cell(0, 6);
for si = 1:numel(subjList)
    subName = subjList{si};
    fprintf('[%s] fNIRS features: %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, si, numel(subjList));

    subDir = fullfile(fnirsRoot, subName);
    inMat = fullfile(subDir, sprintf('%s_task-%s_fnirs_preproc.mat', subName, taskName));
    if exist(inMat, 'file') ~= 2
        tmp = dir(fullfile(subDir, '*_preproc.mat'));
        if isempty(tmp)
            fprintf('[%s] WARNING: fNIRS preproc MAT missing for %s. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
            continue;
        end
        inMat = fullfile(tmp(1).folder, tmp(1).name);
    end

    S = load(inMat, 'dc', 'stim', 'mlActAuto', 'tIncAuto');
    if ~isfield(S, 'dc') || ~isfield(S, 'stim')
        fprintf('[%s] WARNING: Missing dc/stim in %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), inMat);
        continue;
    end

    dc = S.dc;
    stim = S.stim;
    if iscell(stim)
        stim = [stim{:}];
    end

    t = double(dc.time);
    y = double(dc.dataTimeSeries);
    nT = numel(t);
    if size(y, 1) ~= nT
        fprintf('[%s] WARNING: Time length mismatch for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        fprintf('[%s] WARNING: Invalid time step for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    nStim = numel(stim);
    if nStim < max(stimIndexCommon, stimIndexOddball)
        fprintf('[%s] WARNING: stim entries missing for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    stimCommon = stim(stimIndexCommon).data;
    stimOddball = stim(stimIndexOddball).data;
    if isempty(stimCommon) || isempty(stimOddball)
        fprintf('[%s] WARNING: stim common/oddball empty for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    if size(stimCommon, 2) < 3
        stimCommon(:, 3) = 1;
    end
    if size(stimOddball, 2) < 3
        stimOddball(:, 3) = 1;
    end

    % HRF
    tHrf = (0:dt:hrfDurationSec)';
    hrf = (tHrf.^(gammaShape-1)) .* exp(-tHrf/gammaScale) ./ (gamma(gammaShape) * (gammaScale^gammaShape));
    if max(hrf) > 0
        hrf = hrf ./ max(hrf);
    end

    uCommon = zeros(nT, 1);
    for ei = 1:size(stimCommon, 1)
        onset = stimCommon(ei, 1);
        dur = stimCommon(ei, 2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        uCommon(idx) = uCommon(idx) + 1;
    end
    xCommon = conv(uCommon, hrf);
    xCommon = xCommon(1:nT);

    uOdd = zeros(nT, 1);
    for ei = 1:size(stimOddball, 1)
        onset = stimOddball(ei, 1);
        dur = stimOddball(ei, 2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        uOdd(idx) = uOdd(idx) + 1;
    end
    xOdd = conv(uOdd, hrf);
    xOdd = xOdd(1:nT);

    xAll = xCommon + xOdd;

    % Drift regressors
    Xdrift = [];
    switch lower(strtrim(driftModel))
        case 'none'
            % no drift
        case 'linear'
            tTrend = t - mean(t);
            Xdrift = [Xdrift tTrend];
        case 'poly'
            if driftPolyOrder >= 1
                t0 = t - mean(t);
                tNorm = t0 / max(abs(t0));
                for po = 1:driftPolyOrder
                    Xdrift = [Xdrift tNorm.^po]; %#ok<AGROW>
                end
            end
        case 'dct'
            totalDur = (nT - 1) * dt;
            if totalDur > 0 && dctCutoffSec > 0
                K = floor((2 * totalDur) / dctCutoffSec);
                if K >= 1
                    n = (0:nT-1)';
                    dctReg = zeros(nT, K);
                    for k = 1:K
                        dctReg(:, k) = cos(pi * (n + 0.5) * k / nT);
                    end
                    Xdrift = [Xdrift dctReg];
                end
            end
        otherwise
            error('Unknown driftModel: %s', driftModel);
    end

    Xall = [xAll ones(nT,1) Xdrift];
    Xcond = [xCommon xOdd ones(nT,1) Xdrift];

    % Quality masks
    tIncVec = true(nT, 1);
    if applyFnirsTIncAuto && isfield(S, 'tIncAuto') && ~isempty(S.tIncAuto)
        tmp = S.tIncAuto;
        if iscell(tmp)
            tmp = tmp{1};
        end
        tmp = double(tmp(:) > 0);
        if numel(tmp) == nT
            tIncVec = tmp > 0;
        end
    end

    ml = dc.measurementList;
    nMeasAll = numel(ml);
    src = [ml.sourceIndex]';
    det = [ml.detectorIndex]';
    if isprop(ml(1), 'dataTypeLabel')
        chrom = string({ml.dataTypeLabel})';
    else
        dataType = [ml.dataType]';
        chrom = strings(nMeasAll, 1);
        chrom(dataType == 1) = "HbO";
        chrom(dataType == 2) = "HbR";
        chrom(dataType == 3) = "HbT";
        chrom(chrom == "") = "UNK";
    end

    mlAct = true(nMeasAll, 1);
    if isfield(S, 'mlActAuto') && ~isempty(S.mlActAuto)
        tmp = S.mlActAuto;
        if iscell(tmp)
            tmp = tmp{1};
        end
        tmp = double(tmp(:) > 0);
        if numel(tmp) == nMeasAll
            mlAct = tmp > 0;
        end
    end

    for ri = 1:numel(roiDefs)
        roiName = roiDefs(ri).name;
        pairs = roiDefs(ri).pairs;
        roiMask = false(nMeasAll, 1);
        for pi = 1:size(pairs, 1)
            roiMask = roiMask | (src == pairs(pi, 1) & det == pairs(pi, 2));
        end

        for ci = 1:numel(chromList)
            chromName = chromList{ci};
            chromMask = roiMask & mlAct & (chrom == string(chromName));
            nCh = sum(chromMask);
            if nCh < minChannels
                fprintf('[%s] WARNING: No channels for ROI %s %s in %s.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiName, chromName, subName);
                continue;
            end

            yRoi = mean(y(:, chromMask), 2);
            goodMask = tIncVec & isfinite(yRoi) & isfinite(xAll) & isfinite(xCommon) & isfinite(xOdd);
            if sum(goodMask) < size(Xcond, 2) + 2
                fprintf('[%s] WARNING: Too few valid samples for %s %s %s.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, chromName);
                continue;
            end

            bAll = Xall(goodMask, :) \ yRoi(goodMask);
            bCond = Xcond(goodMask, :) \ yRoi(goodMask);
            betaAll = bAll(1);
            betaCommon = bCond(1);
            betaOdd = bCond(2);
            betaDiff = betaOdd - betaCommon;

            fnirsRows(end+1, :) = {subName, roiName, chromName, 'All', betaAll, nCh}; %#ok<SAGROW>
            fnirsRows(end+1, :) = {subName, roiName, chromName, 'Common', betaCommon, nCh}; %#ok<SAGROW>
            fnirsRows(end+1, :) = {subName, roiName, chromName, 'OddballCommon', betaDiff, nCh}; %#ok<SAGROW>
        end
    end
end

if isempty(fnirsRows)
    error('No fNIRS GLM features computed. Check fNIRS paths and settings.');
end

fnirsFeatures = cell2table(fnirsRows, 'VariableNames', ...
    {'SubID', 'ROI', 'Chrom', 'Analysis', 'Beta', 'NChannels'});

%% Alignment summary (oddball trials only)
alignmentSummary = table;
if useAlignmentSummary
    alignRows = cell(0, 8);
    for si = 1:numel(subjList)
        subName = subjList{si};
        alignPath = fullfile(alignmentRoot, subName, ...
            sprintf('%s_task-%s_eeg_fnirs_alignment.csv', subName, taskName));
        if exist(alignPath, 'file') ~= 2
            continue;
        end
        A = readtable(alignPath);
        required = {'eegRejected', 'fnirsBadTime', 'eegP300'};
        if ~all(ismember(required, A.Properties.VariableNames))
            continue;
        end
        eegRejected = A.eegRejected;
        fnirsBad = A.fnirsBadTime;
        eegP300 = A.eegP300;
        oddTotal = height(A);
        eegValid = sum(eegRejected == 0 & ~isnan(eegP300));
        fnirsValid = sum(fnirsBad == 0);
        bothValid = sum(eegRejected == 0 & fnirsBad == 0 & ~isnan(eegP300));
        alignRows(end+1, :) = {subName, oddTotal, eegValid, fnirsValid, bothValid, ...
            sum(eegRejected ~= 0), sum(fnirsBad ~= 0), alignPath}; %#ok<SAGROW>
    end
    if ~isempty(alignRows)
        alignmentSummary = cell2table(alignRows, 'VariableNames', ...
            {'SubID', 'OddballTotal', 'EEGValid', 'FNIRSValid', 'MatchedValid', ...
            'EEGRejected', 'FNIRSBadTime', 'AlignmentFile'});
    end
end

%% Correlation analysis
baseEegFeatures = {'P300_MeanAmp', 'P300_PeakLatencyMs'};
basePowNames = strcat(eegPowerBandNames, 'Power');
baseEegFeatures = [baseEegFeatures, basePowNames];

corrRows = cell(0, 8);
hasTcdf = exist('tcdf', 'file') == 2;
if ~hasTcdf
    fprintf('[%s] WARNING: tcdf not found; p-values will be NaN.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

for ai = 1:numel(analysisNames)
    analysisName = analysisNames{ai};
    suffix = analysisSuffix{ai};

    for ef = 1:numel(baseEegFeatures)
        eegCol = sprintf('%s_%s', baseEegFeatures{ef}, suffix);
        if ~ismember(eegCol, eegFeatures.Properties.VariableNames)
            continue;
        end

        for ri = 1:numel(roiDefs)
            roiName = roiDefs(ri).name;
            for ci = 1:numel(chromList)
                chromName = chromList{ci};
                mask = strcmp(fnirsFeatures.ROI, roiName) & strcmp(fnirsFeatures.Chrom, chromName) & ...
                    strcmp(fnirsFeatures.Analysis, analysisName);
                if ~any(mask)
                    continue;
                end
                subIds = fnirsFeatures.SubID(mask);
                betaVals = fnirsFeatures.Beta(mask);

                [commonSubs, idxFnirs, idxEeg] = intersect(subIds, eegFeatures.SubID, 'stable'); %#ok<ASGLU>
                if isempty(idxEeg)
                    continue;
                end
                x = eegFeatures.(eegCol)(idxEeg);
                y = betaVals(idxFnirs);
                valid = isfinite(x) & isfinite(y);
                n = sum(valid);
                r = NaN;
                p = NaN;
                if n >= minSubjectsForCorr && std(x(valid)) > 0 && std(y(valid)) > 0
                    R = corrcoef(x(valid), y(valid));
                    r = R(1, 2);
                    if hasTcdf
                        tVal = r * sqrt((n - 2) / max(eps, (1 - r^2)));
                        p = 2 * (1 - tcdf(abs(tVal), n - 2));
                    end
                end
                corrRows(end+1, :) = {analysisName, eegCol, roiName, chromName, n, r, p, ...
                    strjoin(commonSubs(valid), ';')}; %#ok<SAGROW>
            end
        end
    end
end

corrTable = cell2table(corrRows, 'VariableNames', ...
    {'Analysis', 'EEGFeature', 'ROI', 'Chrom', 'N', 'R', 'P', 'Subjects'});

%% Scatter plots
if saveScatterPlots && ~isempty(corrTable)
    figDir = fullfile(outRoot, 'figures');
    if exist(figDir, 'dir') ~= 7
        mkdir(figDir);
    end
    plotCount = 0;
    for i = 1:height(corrTable)
        if plotCount >= scatterMaxPlots
            break;
        end
        n = corrTable.N(i);
        if n < minSubjectsForCorr
            continue;
        end
        analysisName = corrTable.Analysis{i};
        eegCol = corrTable.EEGFeature{i};
        roiName = corrTable.ROI{i};
        chromName = corrTable.Chrom{i};

        mask = strcmp(fnirsFeatures.Analysis, analysisName) & ...
            strcmp(fnirsFeatures.ROI, roiName) & strcmp(fnirsFeatures.Chrom, chromName);
        if ~any(mask) || ~ismember(eegCol, eegFeatures.Properties.VariableNames)
            continue;
        end

        subIds = fnirsFeatures.SubID(mask);
        betaVals = fnirsFeatures.Beta(mask);
        [commonSubs, idxFnirs, idxEeg] = intersect(subIds, eegFeatures.SubID, 'stable'); %#ok<ASGLU>
        if isempty(idxEeg)
            continue;
        end
        x = eegFeatures.(eegCol)(idxEeg);
        y = betaVals(idxFnirs);
        valid = isfinite(x) & isfinite(y);
        if sum(valid) < minSubjectsForCorr
            continue;
        end

        plotCount = plotCount + 1;
        fig = figure('Visible', scatterVisible);
        scatter(x(valid), y(valid), 36, 'filled');
        hold on;
        if scatterFitLine && sum(valid) >= 2
            pfit = polyfit(x(valid), y(valid), 1);
            xfit = linspace(min(x(valid)), max(x(valid)), 100);
            yfit = polyval(pfit, xfit);
            plot(xfit, yfit, 'k-', 'LineWidth', 1.2);
        end
        xlabel(eegCol, 'Interpreter', 'none');
        ylabel(sprintf('%s %s %s Beta', analysisName, roiName, chromName), 'Interpreter', 'none');
        title(sprintf('%s | %s | %s | r=%.3f p=%.3f n=%d', ...
            analysisName, eegCol, roiName, corrTable.R(i), corrTable.P(i), n), 'Interpreter', 'none');
        grid on;

        fname = sprintf('scatter_%s_%s_%s_%s.%s', analysisName, eegCol, roiName, chromName, scatterFigureFormat);
        fname = regexprep(fname, '[^A-Za-z0-9._-]', '_');
        fpath = fullfile(figDir, fname);
        if exist('exportgraphics', 'file') == 2
            exportgraphics(fig, fpath, 'Resolution', scatterDpi);
        else
            print(fig, fpath, ['-d' scatterFigureFormat], sprintf('-r%d', scatterDpi));
        end
        close(fig);
    end
end

%% Save outputs
writetable(eegFeatures, fullfile(outRoot, 'STEP1_eeg_features.csv'));
writetable(fnirsFeatures, fullfile(outRoot, 'STEP1_fnirs_features.csv'));
writetable(corrTable, fullfile(outRoot, 'STEP1_group_correlation.csv'));
if ~isempty(alignmentSummary)
    writetable(alignmentSummary, fullfile(outRoot, 'STEP1_alignment_summary.csv'));
end

fprintf('[%s] Done. Outputs in %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), outRoot);
