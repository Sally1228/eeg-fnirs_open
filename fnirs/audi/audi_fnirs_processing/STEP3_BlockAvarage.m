%% STEP2_BlockAvarage.m
% 可选：当成8个blcck，不区分oddball和common sound，按照block-design计算，会丢失一些信息
% 回答的问题：连续听声音+持续做反应”相对休息，是否引出总体血流变化？
% Block-average analysis for auditory oddball task
%
% Assumptions & Notes:
% - stim(1) is common sound (marker=1), stim(2) is oddball sound (marker=2).
% - Block onsets are defined by large gaps in trial onsets (blockGapSec).
% - If event durations are missing, use stimDurationSec as a fallback.
% - Block duration is auto-estimated from blockwise onset span unless set.

clear; clc;

%% ===================== [1] CONFIG =====================
% Optional: Homer3 path
homer3Dir = '';

% Root derived from this script location
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
rootDir = fileparts(thisDir);

% Subjects to process (empty = auto-discover sub-*)
subjList = {};

% Input preproc MAT naming
preprocSuffix = '_task-audi_fnirs_preproc.mat';

% Output folders (per-subject, set inside loop)
outDir = '';
figDir = '';

% Stim indices in preproc MAT
stimIndexCommon = 1;
stimIndexOddball = 2;

% New stim names (labels)
mergedStimName = 'merged12';
blockStimName = 'block_onset';

% Trial/block structure
stimDurationSec = 0.5;
trialISIsec = 2.0;
blockTrials = 23;
expectedTotalTrials = 184;
expectedBlocks = 8;
blockDurationSec = []; % [] = auto from blockwise onset span (median)
blockGapSec = 10; % gap threshold (s) to detect block boundaries

% Epoching for block average
baselineWindowSec = [-5 0];
doBaseline = true;
postBlockSec = 15;

% Plot colors
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue

% ROI definitions (source-detector pairs)
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

fprintf('[INFO] Assumptions: stim(1)=common, stim(2)=oddball, blockGapSec=%.1f\n', blockGapSec);

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

    if nStim < max(stimIndexCommon, stimIndexOddball)
        error('stim has only %d entries; need stim(%d) and stim(%d).', ...
            nStim, stimIndexCommon, stimIndexOddball);
    end

    %% ===================== [5] BUILD MERGED & BLOCK STIM =====================
    stimCommon = stim(stimIndexCommon);
    stimOddball = stim(stimIndexOddball);

    % Extract stim data (ensure Nx3)
    dataCommon = stimCommon.data;
    dataOddball = stimOddball.data;

    if isempty(dataCommon) || isempty(dataOddball)
        error('stim(%d) or stim(%d) has empty data.', stimIndexCommon, stimIndexOddball);
    end

    if size(dataCommon,2) == 1
        dataCommon = [dataCommon, repmat(stimDurationSec, size(dataCommon,1), 1), ones(size(dataCommon,1),1)];
    elseif size(dataCommon,2) == 2
        dataCommon = [dataCommon(:,1), dataCommon(:,2), ones(size(dataCommon,1),1)];
    else
        dataCommon = dataCommon(:,1:3);
    end

    if size(dataOddball,2) == 1
        dataOddball = [dataOddball, repmat(stimDurationSec, size(dataOddball,1), 1), ones(size(dataOddball,1),1)];
    elseif size(dataOddball,2) == 2
        dataOddball = [dataOddball(:,1), dataOddball(:,2), ones(size(dataOddball,1),1)];
    else
        dataOddball = dataOddball(:,1:3);
    end

    dataMerged = [dataCommon; dataOddball];
    [~, sortIdx] = sort(dataMerged(:,1));
    dataMerged = dataMerged(sortIdx,:);

    totalTrials = size(dataMerged, 1);
    if expectedTotalTrials > 0 && totalTrials ~= expectedTotalTrials
        fprintf('[WARN] Total trials=%d (expected %d)\n', totalTrials, expectedTotalTrials);
    end

    onsetsAll = dataMerged(:,1);
    itiAll = diff(onsetsAll); % raw inter-trial gaps for QA
    iti = itiAll; %计算相邻起始时间的差值，得到每个相邻 trial 的时间间隔
    iti = iti(isfinite(iti));% 去掉 NaN/Inf，只保留有效的间隔值。
    if isempty(iti)
        medianITI = trialISIsec;
    else
        medianITI = median(iti);
    end
    % Detect block starts by large gaps in onsets (robust to missing/misaligned trials).
    gapIdx = find(diff(onsetsAll) > blockGapSec);
    blockStartIdx = [1; gapIdx + 1];
    blockOnsets = onsetsAll(blockStartIdx);

    % Block duration: per-block onset span + stim duration, then take median.
    blockEndIdx = [blockStartIdx(2:end) - 1; totalTrials];
    blockDurations = onsetsAll(blockEndIdx) - onsetsAll(blockStartIdx) + stimDurationSec;
    blockDurations = blockDurations(isfinite(blockDurations) & blockDurations > 0);

    blockDurationSecLocal = blockDurationSec;
    if isempty(blockDurationSecLocal)
        if ~isempty(blockDurations)
            blockDurationSecLocal = median(blockDurations);
        else
            blockDurationSecLocal = blockTrials * medianITI;
        end
    end
    blockDurationSecUsed = blockDurationSecLocal;

    if expectedBlocks > 0 && numel(blockOnsets) ~= expectedBlocks
        fprintf('[WARN] Blocks=%d (expected %d)\n', numel(blockOnsets), expectedBlocks);
    end

    dataBlock = [blockOnsets, repmat(blockDurationSecUsed, numel(blockOnsets), 1), ones(numel(blockOnsets),1)];

    % Build new stim entries
    stimMerged = stimCommon;
    stimMerged.name = mergedStimName;
    stimMerged.data = dataMerged;

    stimBlock = stimCommon;
    stimBlock.name = blockStimName;
    stimBlock.data = dataBlock;

    % Insert / replace by name
    mergedIdx = find(strcmp(stimNames, mergedStimName), 1);
    if isempty(mergedIdx)
        stim = [stim(:); stimMerged];
    else
        stim(mergedIdx) = stimMerged;
    end

    blockIdx = find(strcmp(stimNames, blockStimName), 1);
    if isempty(blockIdx)
        stim = [stim(:); stimBlock];
    else
        stim(blockIdx) = stimBlock;
    end

    fprintf('[INFO] Added/updated stim: %s (trials=%d), %s (blocks=%d)\n', ...
        mergedStimName, totalTrials, blockStimName, numel(blockOnsets));

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

    epochStart = baselineWindowSec(1);
    epochEnd = blockDurationSecUsed + postBlockSec;
    epochTime = (epochStart:dt:epochEnd)';

    baseMask = epochTime >= baselineWindowSec(1) & epochTime <= baselineWindowSec(2);
    blockMask = epochTime >= 0 & epochTime <= blockDurationSecUsed;

    blockCount = numel(blockOnsets);
    results = struct();
    results.epochTime = epochTime;
    results.blockOnsets = blockOnsets;
    results.blockDurationSec = blockDurationSecUsed;
    results.blockDurations = blockDurations;
    results.blockStartIdx = blockStartIdx;
    results.blockGapSec = blockGapSec;
    results.roi = struct();

    rowCell = {};
    for ri = 1:nRoi
        hboMat = nan(numel(epochTime), blockCount);
        hbrMat = nan(numel(epochTime), blockCount);
        hboMeans = nan(blockCount,1);
        hbrMeans = nan(blockCount,1);

        for bi = 1:blockCount
            tQuery = blockOnsets(bi) + epochTime;
            hboTmp = interp1(t, roiHbo{ri}, tQuery, 'linear', NaN);
            hbrTmp = interp1(t, roiHbr{ri}, tQuery, 'linear', NaN);

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

        results.roi(ri).name = roiDefs(ri).name;
        results.roi(ri).hboAvg = mean(hboMat, 2, 'omitnan');
        results.roi(ri).hbrAvg = mean(hbrMat, 2, 'omitnan');
        results.roi(ri).hboSem = hboSem;
        results.roi(ri).hbrSem = hbrSem;
        results.roi(ri).hboCi = hboCi;
        results.roi(ri).hbrCi = hbrCi;
        results.roi(ri).hboBlockMeans = hboMeans;
        results.roi(ri).hbrBlockMeans = hbrMeans;
        results.roi(ri).nHboChannels = roiNhbo(ri);
        results.roi(ri).nHbrChannels = roiNhbr(ri);

        for bi = 1:blockCount
            rowCell(end+1,:) = {subName, roiDefs(ri).name, 'HbO', bi, hboMeans(bi), roiNhbo(ri)}; %#ok<AGROW>
            rowCell(end+1,:) = {subName, roiDefs(ri).name, 'HbR', bi, hbrMeans(bi), roiNhbr(ri)}; %#ok<AGROW>
        end
    end

    %% ===================== [8] SAVE OUTPUTS =====================
    outMat = fullfile(outDir, sprintf('%s_task-audi_blockavg.mat', subName));
    save(outMat, 'results', 'stim', 'roiDefs', 'blockOnsets', 'blockDurationSecUsed', '-v7.3');

    outCsv = fullfile(outDir, sprintf('%s_task-audi_blockavg_table.csv', subName));
    T = cell2table(rowCell, 'VariableNames', ...
        {'subjId','roi','chromophore','blockIndex','blockMean','nChannels'});
    writetable(T, outCsv);

    %% ===================== [9] FIGURES =====================
    % QA: block boundary detection
    fig0 = figure('Visible','off');
    subplot(2,1,1);
    plot(1:totalTrials, onsetsAll, 'k.-'); hold on;
    plot(blockStartIdx, blockOnsets, 'ro', 'MarkerFaceColor','r');
    ylabel('Onset (s)');
    title(sprintf('%s block boundary QA', subName), 'Interpreter','none');
    legend({'trial onsets','block starts'}, 'Location','best');
    grid on;

    subplot(2,1,2);
    if ~isempty(itiAll)
        gapX = 1:numel(itiAll);
        plot(gapX, itiAll, 'k.-'); hold on;
        plot([gapX(1) gapX(end)], [blockGapSec blockGapSec], '--r');
        if ~isempty(gapIdx)
            plot(gapIdx, itiAll(gapIdx), 'ro', 'MarkerFaceColor','r');
        end
        xlabel('Trial index (gap from i to i+1)');
        ylabel('Inter-trial gap (s)');
        legend({'gap','threshold','gap>threshold'}, 'Location','best');
        grid on;
    else
        text(0.5, 0.5, 'Not enough trials for gap plot', 'HorizontalAlignment','center');
        axis off;
    end
    figPath0 = fullfile(figDir, sprintf('%s_task-audi_step-blockavg_qc_blockboundaries.png', subName));
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
    figPath1 = fullfile(figDir, sprintf('%s_task-audi_step-blockavg_qc_channels.png', subName));
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
    figPath2 = fullfile(figDir, sprintf('%s_task-audi_step-blockavg_signal.png', subName));
    saveas(fig2, figPath2);
    close(fig2);

    % Result plot: block average time course
    fig3 = figure('Visible','off');
    for ri = 1:nRoi
        subplot(nRoi,1,ri);
        x = results.epochTime;
        hboAvg = results.roi(ri).hboAvg;
        hbrAvg = results.roi(ri).hbrAvg;
        hboCi = results.roi(ri).hboCi;
        hbrCi = results.roi(ri).hbrCi;
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
        ylabel(roiDefs(ri).name, 'Interpreter','none');
        if ri == 1
            title(sprintf('%s block average (baseline-corrected)', subName), 'Interpreter','none');
        end
        if ri == nRoi
            xlabel('Time (s)');
        end
    end
    figPath3 = fullfile(figDir, sprintf('%s_task-audi_step-blockavg_result.png', subName));
    saveas(fig3, figPath3);
    close(fig3);

    % Trend plot: block index vs block mean (drift/fatigue check)
    fig4 = figure('Visible','off');
    for ri = 1:nRoi
        subplot(nRoi,1,ri);
        x = 1:blockCount;
        plot(x, results.roi(ri).hboBlockMeans, '-o', 'Color', colorHbO, 'LineWidth', 1); hold on;
        plot(x, results.roi(ri).hbrBlockMeans, '-o', 'Color', colorHbR, 'LineWidth', 1);
        ylabel(roiDefs(ri).name, 'Interpreter','none');
        if ri == 1
            title(sprintf('%s block mean trend', subName), 'Interpreter','none');
            legend({'HbO','HbR'}, 'Location','best');
        end
        if ri == nRoi
            xlabel('Block index');
        end
        grid on;
    end
    figPath4 = fullfile(figDir, sprintf('%s_task-audi_step-blockavg_trend.png', subName));
    saveas(fig4, figPath4);
    close(fig4);

    fprintf('[DONE] %s | blocks=%d | out=%s\n', subName, blockCount, outMat);

    if processOnlyFirst
        break;
    end
end
