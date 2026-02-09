%% STEP3_BlockAvarage.m
% 可选：8个block，其中4个左手，4个右手
% Assumptions & Notes:
% - stim(1)=left hand, stim(2)=right hand.
% - Left uses RightCentral ROI; right uses LeftCentral ROI.
% - Block boundaries are detected by blockGapSec.
clear; clc;

%% ===================== [1] CONFIG =====================
% Optional: Homer3 path
homer3Dir = '';

% Root derived from this script location
% thisDir ='/Users/zhaifeifei/Desktop/eeg_fnirs/fnirs/motor/motor_fnirs_processing'
thisDir = fileparts(mfilename('fullpath'));

rootDir = fileparts(thisDir);

% Subjects to process (empty = auto-discover sub-*)
subjList = {};

% Input preproc MAT naming
preprocSuffix = '_task-motor_fnirs_preproc.mat';

% Output folders (per-subject, set inside loop)
outDir = '';
figDir = '';

% Stim indices in preproc MAT (1=left, 2=right)
stimIndexLeft = 1;
stimIndexRight = 2;

% No merged/block stim output; keep original stim(1) and stim(2) only.
% Contralateral ROI mapping
contraRoiForLeft = 'RightCentral';
contraRoiForRight = 'LeftCentral';

% Trial/block structure (block design: 30s task + 17s rest)
stimDurationSec = 30;
trialISIsec = 47;
blockTrials = 1;
expectedTotalTrials = 8;
expectedBlocksPerCond = 4;
blockDurationSec = []; % [] = auto from blockwise onset span (median)
blockGapSec = 20; % gap threshold (s) to detect block boundaries

% Epoching for block average
baselineWindowSec = [-5 0];
doBaseline = true;
postBlockSec = 10;

% Plot colors
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue

% ROI definitions (source-detector pairs)
%注意检查一下左右是不是正确？可以再做一下对侧

roiDefs = struct();
roiDefs(1).name = 'LeftCentral';
roiDefs(1).pairs = [ ...
    14 10; 10 10; 10 14; 14 14; 14 18; 18 18; 18 14];
roiDefs(2).name = 'RightCentral';
roiDefs(2).pairs = [ ...
    11 11; 15 11; 15 15; 11 15; 19 15; 19 19; 15 19];

% Debug
processOnlyFirst = false;

%% ===================== [2] ADD HOMER3 TO PATH (OPTIONAL) =====================
if ~isempty(homer3Dir)
    if exist(homer3Dir, 'dir') ~= 7
        error('homer3Dir does not exist: %s', homer3Dir);
    end
    addpath(genpath(homer3Dir));
end

%% ===================== [3] DISCOVER INPUT FILES =====================
if isempty(subjList)
    subDirs = dir(fullfile(rootDir, 'sub-*'));
    if isempty(subDirs)
        error('No sub-* folders found under: %s', rootDir);
    end
    subjList = {subDirs.name};
end

fprintf('[INFO] Assumptions: stim(1)=left->%s, stim(2)=right->%s, blockGapSec=%.1f\n', ...
    contraRoiForLeft, contraRoiForRight, blockGapSec);

%% ===================== [4] MAIN LOOP =====================
for si = 1:numel(subjList)
    subName = subjList{si};
    subPath = fullfile(rootDir, subName);
    inMat = fullfile(subPath, [subName preprocSuffix]);

    if exist(inMat, 'file') ~= 2
        % fallback: pick any *_preproc.mat in subject folder
        tmp = dir(fullfile(subPath, '*_preproc.mat'));
        if isempty(tmp)
            fprintf('[SKIP] %s: no preproc MAT found.\n', subName);
            continue;
        end
        inMat = fullfile(tmp(1).folder, tmp(1).name);
    end

    fprintf('\n==================== SUBJECT: %s ====================\n', subName);
    fprintf('Loading: %s\n', inMat);

    outDir = fullfile(subPath, 'blockavg');
    figDir = outDir;
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    if ~exist(figDir, 'dir'); mkdir(figDir); end

    S = load(inMat, 'dc', 'stim', 'qc');
    try
        Sml = load(inMat, 'mlActAuto');
        S.mlActAuto = Sml.mlActAuto;
    catch
    end
    try
        Sdod = load(inMat, 'dod');
        S.dod = Sdod.dod;
    catch
    end
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

    % --- data & time ---
    if ~isprop(dc, 'dataTimeSeries') && ~isfield(dc, 'dataTimeSeries')
        error('dc has no dataTimeSeries field.');
    end
    if ~isprop(dc, 'time') && ~isfield(dc, 'time')
        error('dc has no time field.');
    end
    t = double(dc.time);
    y = double(dc.dataTimeSeries);
    if size(y,1) ~= numel(t)
        error('Time length (%d) does not match data rows (%d).', numel(t), size(y,1));
    end

    % --- stim summary ---
    nStim = numel(stim);
    fprintf('[INFO] Stim summary (%d entries):\n', nStim);
    stimNames = cell(nStim,1);
    for k = 1:nStim
        if isprop(stim(k), 'name')
            stimNames{k} = char(stim(k).name);
        elseif isfield(stim(k), 'name')
            stimNames{k} = char(stim(k).name);
        else
            stimNames{k} = sprintf('stim(%d)', k);
        end

        if isprop(stim(k), 'data')
            stimData = stim(k).data;
        else
            stimData = stim(k).data;
        end
        if isempty(stimData)
            nEv = 0;
        else
            nEv = size(stimData, 1);
        end
        fprintf('  - %s: %d events\n', stimNames{k}, nEv);
    end

    if nStim < max(stimIndexLeft, stimIndexRight)
        error('stim has only %d entries; need stim(%d) and stim(%d).', ...
            nStim, stimIndexLeft, stimIndexRight);
    end

    %% ===================== [5] BUILD CONDITIONS =====================
    stimLeft = stim(stimIndexLeft);
    stimRight = stim(stimIndexRight);

    % Extract stim data (ensure Nx3)
    dataLeft = stimLeft.data;
    dataRight = stimRight.data;

    if isempty(dataLeft) || isempty(dataRight)
        error('stim(%d) or stim(%d) has empty data.', stimIndexLeft, stimIndexRight);
    end

    if size(dataLeft,2) == 1
        dataLeft = [dataLeft, repmat(stimDurationSec, size(dataLeft,1), 1), ones(size(dataLeft,1),1)];
    elseif size(dataLeft,2) == 2
        dataLeft = [dataLeft(:,1), dataLeft(:,2), ones(size(dataLeft,1),1)];
    else
        dataLeft = dataLeft(:,1:3);
    end

    if size(dataRight,2) == 1
        dataRight = [dataRight, repmat(stimDurationSec, size(dataRight,1), 1), ones(size(dataRight,1),1)];
    elseif size(dataRight,2) == 2
        dataRight = [dataRight(:,1), dataRight(:,2), ones(size(dataRight,1),1)];
    else
        dataRight = dataRight(:,1:3);
    end

    dataLeft = sortrows(dataLeft, 1);
    dataRight = sortrows(dataRight, 1);

    totalTrials = size(dataLeft, 1) + size(dataRight, 1);
    if expectedTotalTrials > 0 && totalTrials ~= expectedTotalTrials
        fprintf('[WARN] Total trials=%d (expected %d)\n', totalTrials, expectedTotalTrials);
    end

    condList = struct();
    condList(1).name = 'left';
    condList(1).data = dataLeft;
    condList(1).contraRoi = contraRoiForLeft;
    condList(2).name = 'right';
    condList(2).data = dataRight;
    condList(2).contraRoi = contraRoiForRight;

    nCond = numel(condList);
    for ci = 1:nCond
        dataCond = condList(ci).data;
        onsets = dataCond(:,1);
        onsets = onsets(:);

        itiAll = diff(onsets); % raw inter-trial gaps for QA
        iti = itiAll;
        iti = iti(isfinite(iti)); % 去掉 NaN/Inf，只保留有效的间隔值
        if isempty(iti)
            medianITI = trialISIsec;
        else
            medianITI = median(iti);
        end

        % Detect block starts by large gaps in onsets (robust to missing/misaligned trials).
        gapIdx = find(diff(onsets) > blockGapSec);
        blockStartIdx = [1; gapIdx + 1];
        blockOnsets = onsets(blockStartIdx);

        % Block duration: per-block onset span + stim duration, then take median.
        blockEndIdx = [blockStartIdx(2:end) - 1; numel(onsets)];
        blockDurations = onsets(blockEndIdx) - onsets(blockStartIdx) + stimDurationSec;
        blockDurations = blockDurations(isfinite(blockDurations) & blockDurations > 0);

        blockDurationSecLocal = blockDurationSec;
        if isempty(blockDurationSecLocal)
            if ~isempty(blockDurations)
                blockDurationSecLocal = median(blockDurations);
            else
                blockDurationSecLocal = blockTrials * medianITI;
            end
        end

        condList(ci).onsets = onsets;
        condList(ci).itiAll = itiAll;
        condList(ci).gapIdx = gapIdx;
        condList(ci).blockStartIdx = blockStartIdx;
        condList(ci).blockOnsets = blockOnsets;
        condList(ci).blockDurations = blockDurations;
        condList(ci).blockDurationSecUsed = blockDurationSecLocal;
        condList(ci).blockCount = numel(blockOnsets);

        if expectedBlocksPerCond > 0 && numel(blockOnsets) ~= expectedBlocksPerCond
            fprintf('[WARN] %s blocks=%d (expected %d)\n', ...
                condList(ci).name, numel(blockOnsets), expectedBlocksPerCond);
        end
        if condList(ci).blockCount == 0
            error('No blocks found for condition %s.', condList(ci).name);
        end
    end

    fprintf('[INFO] Using stim(1) and stim(2) separately; no merged block average.\n');

    %% ===================== [6] ROI CHANNEL SELECTION =====================
    ml = dc.measurementList;
    nMeas = numel(ml);
    if nMeas == 0
        error('dc.measurementList is empty.');
    end

    if isprop(ml(1), 'sourceIndex')
        src = [ml.sourceIndex]';
        det = [ml.detectorIndex]';
    else
        error('measurementList missing sourceIndex/detectorIndex.');
    end

    hboMaskGlobal = [];
    hbrMaskGlobal = [];
    if isprop(ml(1), 'dataTypeLabel')
        dtLabel = {ml.dataTypeLabel}';
        dtLabel = string(dtLabel);
        hboMaskGlobal = dtLabel == "HbO";
        hbrMaskGlobal = dtLabel == "HbR";
    elseif isprop(ml(1), 'dataType')
        dataType = [ml.dataType]';
        hboMaskGlobal = dataType == 1;
        hbrMaskGlobal = dataType == 2;
    else
        error('measurementList missing dataType/dataTypeLabel (HbO/HbR labels).');
    end

    % Apply QC mask (mlActAuto) to exclude bad channels from ROI averages.
    qcMaskGlobal = true(nMeas,1);
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
        dodMl = [];
        if isfield(S, 'dod') && ~isempty(S.dod) && isprop(S.dod, 'measurementList')
            dodMl = S.dod.measurementList;
        end

        if ~isempty(dodMl) && numel(dodMl) == numel(mlActVec) && ...
                isprop(dodMl(1), 'sourceIndex') && isprop(dodMl(1), 'detectorIndex')
            srcDod = [dodMl.sourceIndex]';
            detDod = [dodMl.detectorIndex]';
            pairKeyDod = arrayfun(@(a,b) sprintf('%d_%d', a, b), srcDod, detDod, 'UniformOutput', false);
            [pairKeyDodUniq, ~, pairIdxDod] = unique(pairKeyDod, 'stable');
            pairActive = accumarray(pairIdxDod, mlActVec, [], @all);

            pairKeyDc = arrayfun(@(a,b) sprintf('%d_%d', a, b), src, det, 'UniformOutput', false);
            [pairHit, pairLoc] = ismember(pairKeyDc, pairKeyDodUniq);
            qcMaskGlobal = false(nMeas,1);
            qcMaskGlobal(pairHit) = pairActive(pairLoc(pairHit));
            qcNote = 'mlActAuto mapped by src-det';
            if any(~pairHit)
                fprintf('[WARN] QC mapping: %d dc channels missing from dod ML.\n', sum(~pairHit));
            end
        elseif numel(mlActVec) == nMeas
            qcMaskGlobal = mlActVec;
            qcNote = 'mlActAuto applied directly to dc';
        else
            fprintf('[WARN] mlActAuto size mismatch; ROI will use all channels.\n');
        end
    else
        fprintf('[WARN] mlActAuto not found; ROI will use all channels.\n');
    end
    fprintf('[INFO] QC mask for ROI: keep %d/%d channels (%s)\n', sum(qcMaskGlobal), nMeas, qcNote);

    roiNames = {roiDefs.name};
    nRoi = numel(roiDefs);
    roiHbo = cell(nRoi,1);
    roiHbr = cell(nRoi,1);
    roiNhbo = zeros(nRoi,1);
    roiNhbr = zeros(nRoi,1);

    for ri = 1:nRoi
        pairs = roiDefs(ri).pairs;
        pairMask = false(nMeas,1);
        for pi = 1:size(pairs,1)
            pairMask = pairMask | (src==pairs(pi,1) & det==pairs(pi,2));
        end

        hboIdx = find(pairMask & hboMaskGlobal & qcMaskGlobal);
        hbrIdx = find(pairMask & hbrMaskGlobal & qcMaskGlobal);
        roiNhbo(ri) = numel(hboIdx);
        roiNhbr(ri) = numel(hbrIdx);

        if isempty(hboIdx)
            roiHbo{ri} = nan(numel(t),1);
            fprintf('[WARN] ROI %s: no HbO channels matched.\n', roiDefs(ri).name);
        else
            roiHbo{ri} = mean(y(:,hboIdx), 2, 'omitnan');%求roi平均值
        end

        if isempty(hbrIdx)
            roiHbr{ri} = nan(numel(t),1);
            fprintf('[WARN] ROI %s: no HbR channels matched.\n', roiDefs(ri).name);
        else
            roiHbr{ri} = mean(y(:,hbrIdx), 2, 'omitnan');%求roi平均值
        end
    end

    %% ===================== [7] BLOCK AVERAGE (EPOCH) =====================
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('Invalid time vector; cannot compute epoch sampling interval.');
    end

    results = struct();
    results.conditions = struct();

    rowCell = {};
    for ci = 1:nCond
        condName = condList(ci).name;
        roiIdx = find(strcmp(roiNames, condList(ci).contraRoi), 1);
        if isempty(roiIdx)
            error('Contra ROI not found for %s: %s', condName, condList(ci).contraRoi);
        end

        epochStart = baselineWindowSec(1);
        epochEnd = condList(ci).blockDurationSecUsed + postBlockSec;
        epochTime = (epochStart:dt:epochEnd)';

        baseMask = epochTime >= baselineWindowSec(1) & epochTime <= baselineWindowSec(2);
        blockMask = epochTime >= 0 & epochTime <= condList(ci).blockDurationSecUsed;

        blockCount = condList(ci).blockCount;
        hboMat = nan(numel(epochTime), blockCount);
        hbrMat = nan(numel(epochTime), blockCount);
        hboMeans = nan(blockCount,1);
        hbrMeans = nan(blockCount,1);

        for bi = 1:blockCount
            tQuery = condList(ci).blockOnsets(bi) + epochTime;
            hboTmp = interp1(t, roiHbo{roiIdx}, tQuery, 'linear', NaN);
            hbrTmp = interp1(t, roiHbr{roiIdx}, tQuery, 'linear', NaN);

            if doBaseline
                hboBase = mean(hboTmp(baseMask), 'omitnan');
                hbrBase = mean(hbrTmp(baseMask), 'omitnan');
                hboTmp = hboTmp - hboBase;
                hbrTmp = hbrTmp - hbrBase;
            end

            hboMat(:,bi) = hboTmp;
            hbrMat(:,bi) = hbrTmp;
            hboMeans(bi) = mean(hboTmp(blockMask), 'omitnan');
            hbrMeans(bi) = mean(hbrTmp(blockMask), 'omitnan');
        end

        hboStd = std(hboMat, 0, 2, 'omitnan');
        hbrStd = std(hbrMat, 0, 2, 'omitnan');
        hboN = sum(isfinite(hboMat), 2);
        hbrN = sum(isfinite(hbrMat), 2);
        hboDen = sqrt(hboN); hboDen(hboDen == 0) = 1;
        hbrDen = sqrt(hbrN); hbrDen(hbrDen == 0) = 1;
        hboSem = hboStd ./ hboDen;
        hbrSem = hbrStd ./ hbrDen;
        hboCi = 1.96 * hboSem;
        hbrCi = 1.96 * hbrSem;

        results.conditions(ci).name = condName;
        results.conditions(ci).roiName = roiDefs(roiIdx).name;
        results.conditions(ci).epochTime = epochTime;
        results.conditions(ci).blockOnsets = condList(ci).blockOnsets;
        results.conditions(ci).blockDurationSec = condList(ci).blockDurationSecUsed;
        results.conditions(ci).blockDurations = condList(ci).blockDurations;
        results.conditions(ci).blockStartIdx = condList(ci).blockStartIdx;
        results.conditions(ci).blockGapSec = blockGapSec;
        results.conditions(ci).blockCount = blockCount;
        results.conditions(ci).hboAvg = mean(hboMat, 2, 'omitnan');
        results.conditions(ci).hbrAvg = mean(hbrMat, 2, 'omitnan');
        results.conditions(ci).hboSem = hboSem;
        results.conditions(ci).hbrSem = hbrSem;
        results.conditions(ci).hboCi = hboCi;
        results.conditions(ci).hbrCi = hbrCi;
        results.conditions(ci).hboBlockMeans = hboMeans;
        results.conditions(ci).hbrBlockMeans = hbrMeans;
        results.conditions(ci).nHboChannels = roiNhbo(roiIdx);
        results.conditions(ci).nHbrChannels = roiNhbr(roiIdx);

        for bi = 1:blockCount
            rowCell(end+1,:) = {subName, condName, roiDefs(roiIdx).name, 'HbO', bi, hboMeans(bi), roiNhbo(roiIdx)}; %#ok<AGROW>
            rowCell(end+1,:) = {subName, condName, roiDefs(roiIdx).name, 'HbR', bi, hbrMeans(bi), roiNhbr(roiIdx)}; %#ok<AGROW>
        end
    end

    %% ===================== [8] SAVE OUTPUTS =====================
    outMat = fullfile(outDir, sprintf('%s_task-motor_blockavg.mat', subName));
    save(outMat, 'results', 'stim', 'roiDefs', '-v7.3');

    outCsv = fullfile(outDir, sprintf('%s_task-motor_blockavg_table.csv', subName));
    T = cell2table(rowCell, 'VariableNames', ...
        {'subjId','condition','roi','chromophore','blockIndex','blockMean','nChannels'});
    writetable(T, outCsv);

    %% ===================== [9] FIGURES =====================
    % QA: block boundary detection
    fig0 = figure('Visible','off');
    for ci = 1:nCond
        subplot(nCond,2,(ci-1)*2+1);
        plot(1:numel(condList(ci).onsets), condList(ci).onsets, 'k.-'); hold on;
        plot(condList(ci).blockStartIdx, condList(ci).blockOnsets, 'ro', 'MarkerFaceColor','r');
        ylabel(sprintf('%s onset (s)', condList(ci).name), 'Interpreter','none');
        title(sprintf('%s block boundary QA (%s)', subName, condList(ci).name), 'Interpreter','none');
        grid on;

        subplot(nCond,2,(ci-1)*2+2);
        if ~isempty(condList(ci).itiAll)
            gapX = 1:numel(condList(ci).itiAll);
            plot(gapX, condList(ci).itiAll, 'k.-'); hold on;
            plot([gapX(1) gapX(end)], [blockGapSec blockGapSec], '--r');
            if ~isempty(condList(ci).gapIdx)
                plot(condList(ci).gapIdx, condList(ci).itiAll(condList(ci).gapIdx), 'ro', 'MarkerFaceColor','r');
            end
            xlabel('Trial index (gap from i to i+1)');
            ylabel('Inter-trial gap (s)');
            legend({'gap','threshold','gap>threshold'}, 'Location','best');
            grid on;
        else
            text(0.5, 0.5, 'Not enough trials for gap plot', 'HorizontalAlignment','center');
            axis off;
        end
    end
    figPath0 = fullfile(figDir, sprintf('%s_task-motor_step-blockavg_qc_blockboundaries.png', subName));
    saveas(fig0, figPath0);
    close(fig0);

    % QC: channel counts per ROI
    fig1 = figure('Visible','off');
    b = bar(categorical(roiNames), [roiNhbo roiNhbr]);
    if numel(b) >= 1
        b(1).FaceColor = colorHbO;
    end
    if numel(b) >= 2
        b(2).FaceColor = colorHbR;
    end
    ylabel('Channels');
    title(sprintf('%s ROI channel counts', subName), 'Interpreter','none');
    legend({'HbO','HbR'}, 'Location','best');
    figPath1 = fullfile(figDir, sprintf('%s_task-motor_step-blockavg_qc_channels.png', subName));
    saveas(fig1, figPath1);
    close(fig1);

    % Signal plot: ROI mean time series
    fig2 = figure('Visible','off');
    for ri = 1:nRoi
        subplot(nRoi,1,ri);
        plot(t, roiHbo{ri}, 'Color', colorHbO); hold on;
        plot(t, roiHbr{ri}, 'Color', colorHbR);
        ylabel(roiDefs(ri).name, 'Interpreter','none');
        if ri == 1
            title(sprintf('%s ROI mean time series', subName), 'Interpreter','none');
        end
        if ri == nRoi
            xlabel('Time (s)');
        end
    end
    figPath2 = fullfile(figDir, sprintf('%s_task-motor_step-blockavg_signal.png', subName));
    saveas(fig2, figPath2);
    close(fig2);

    % Result plot: block average time course
    fig3 = figure('Visible','off');
    for ci = 1:nCond
        subplot(nCond,1,ci);
        x = results.conditions(ci).epochTime;
        hboAvg = results.conditions(ci).hboAvg;
        hbrAvg = results.conditions(ci).hbrAvg;
        hboCi = results.conditions(ci).hboCi;
        hbrCi = results.conditions(ci).hbrCi;
        hold on;
        if any(isfinite(hboCi))
            fill([x; flipud(x)], [hboAvg + hboCi; flipud(hboAvg - hboCi)], ...
                colorHbO, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        if any(isfinite(hbrCi))
            fill([x; flipud(x)], [hbrAvg + hbrCi; flipud(hbrAvg - hbrCi)], ...
                colorHbR, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        plot(x, hboAvg, 'Color', colorHbO, 'LineWidth', 1.5);
        plot(x, hbrAvg, 'Color', colorHbR, 'LineWidth', 1.5);
        xline(0, '--k');
        ylabel(sprintf('%s | %s', results.conditions(ci).name, results.conditions(ci).roiName), 'Interpreter','none');
        if ci == 1
            title(sprintf('%s block average (baseline-corrected)', subName), 'Interpreter','none');
        end
        if ci == nCond
            xlabel('Time (s)');
        end
    end
    figPath3 = fullfile(figDir, sprintf('%s_task-motor_step-blockavg_result.png', subName));
    saveas(fig3, figPath3);
    close(fig3);

    % Trend plot: block index vs block mean (drift/fatigue check)
    fig4 = figure('Visible','off');
    for ci = 1:nCond
        subplot(nCond,1,ci);
        x = 1:results.conditions(ci).blockCount;
        plot(x, results.conditions(ci).hboBlockMeans, '-o', 'Color', colorHbO, 'LineWidth', 1); hold on;
        plot(x, results.conditions(ci).hbrBlockMeans, '-o', 'Color', colorHbR, 'LineWidth', 1);
        ylabel(sprintf('%s | %s', results.conditions(ci).name, results.conditions(ci).roiName), 'Interpreter','none');
        if ci == 1
            title(sprintf('%s block mean trend', subName), 'Interpreter','none');
            legend({'HbO','HbR'}, 'Location','best');
        end
        if ci == nCond
            xlabel('Block index');
        end
        grid on;
    end
    figPath4 = fullfile(figDir, sprintf('%s_task-motor_step-blockavg_trend.png', subName));
    saveas(fig4, figPath4);
    close(fig4);

    fprintf('[DONE] %s | left=%d right=%d | out=%s\n', subName, ...
        condList(1).blockCount, condList(2).blockCount, outMat);

    if processOnlyFirst
        break;
    end
end
