%% STEP3_BlockAvarage.m
% 语义词语流畅性任务：6个block，marker=1
% Assumptions & Notes:
% - stim(1)=task blocks (marker=1).
% - Task duration=60s, rest=17s between blocks.
% - Block order: Vegetable, Animal, Flower, Instrument, Fruit, Appliance.
% - Block boundaries are detected by blockGapSec.
clear; clc;

%% ===================== [1] CONFIG =====================
% Optional: Homer3 path
homer3Dir = '';

% Root derived from this script location
% thisDir ='/Users/zhaifeifei/Desktop/eeg_fnirs/fnirs/cate/cate_fnirs_processing'
thisDir = fileparts(mfilename('fullpath'));

rootDir = fileparts(thisDir);

% Subjects to process (empty = auto-discover sub-*)
subjList = {};

% Input preproc MAT naming
preprocSuffix = '_task-cate_fnirs_preproc.mat';

% Output folders (per-subject, set inside loop)
outDir = '';
figDir = '';

% Stim index in preproc MAT (1=task marker)
stimIndexTask = 1;

% Trial/block structure (block design: 60s task + 17s rest)
stimDurationSec = 60;
trialISIsec = 77;
blockTrials = 1;
expectedTotalTrials = 6;
expectedBlocksPerCond = 6;
blockDurationSec = []; % [] = auto from blockwise onset span (median)
blockGapSec = 30; % gap threshold (s) to detect block boundaries

% Epoching for block average
baselineWindowSec = [-5 0];
doBaseline = true;
postBlockSec = 17;

% Plot colors
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue

% Block labels (order assumed)
blockLabels = {'Vegetable','Animal','Flower','Instrument','Fruit','Appliance'};

% ROI definitions (source-detector pairs)
roiDefs = struct();
roiDefs(1).name = 'LeftFrontal';
roiDefs(1).pairs = [ ...
    2 2; 5 2; 5 5; 2 5; 5 8; 7 8];
roiDefs(2).name = 'LeftTemporal';
roiDefs(2).pairs = [ ...
    9 9; 9 14; 13 13; 13 9; 17 13; 17 17; 13 17];
roiDefs(3).name = 'RightFrontal';
roiDefs(3).pairs = [ ...
    3 3; 3 6; 6 6; 6 3; 8 6; 8 7];
roiDefs(4).name = 'RightTemporal';
roiDefs(4).pairs = [ ...
    12 16; 12 12; 16 12; 16 16; 20 16; 20 20; 16 20];

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

fprintf('[INFO] Assumptions: stim(1)=marker=1 task blocks, blockGapSec=%.1f\n', blockGapSec);
fprintf('[INFO] Assumptions: task=60s, rest=17s, order=Vegetable->Animal->Flower->Instrument->Fruit->Appliance\n');

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

    if nStim < stimIndexTask
        error('stim has only %d entries; need stim(%d).', nStim, stimIndexTask);
    end

    %% ===================== [5] BUILD CONDITIONS =====================
    stimTask = stim(stimIndexTask);

    % Extract stim data (ensure Nx3)
    dataTask = stimTask.data;

    if isempty(dataTask)
        error('stim(%d) has empty data.', stimIndexTask);
    end

    if size(dataTask,2) == 1
        dataTask = [dataTask, repmat(stimDurationSec, size(dataTask,1), 1), ones(size(dataTask,1),1)];
    elseif size(dataTask,2) == 2
        dataTask = [dataTask(:,1), dataTask(:,2), ones(size(dataTask,1),1)];
    else
        dataTask = dataTask(:,1:3);
    end

    dataTask = sortrows(dataTask, 1);

    totalTrials = size(dataTask, 1);
    if expectedTotalTrials > 0 && totalTrials ~= expectedTotalTrials
        fprintf('[WARN] Total trials=%d (expected %d)\n', totalTrials, expectedTotalTrials);
    end

    condList = struct();
    condList(1).name = 'task';
    condList(1).data = dataTask;
    condList(1).blockLabels = blockLabels;

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

    fprintf('[INFO] Using stim(1) only; no merged block average.\n');

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
    results.task = struct();
    results.rois = struct();

    rowCell = {};
    for ci = 1:nCond
        condName = condList(ci).name;

        epochStart = baselineWindowSec(1);
        epochEnd = condList(ci).blockDurationSecUsed + postBlockSec;
        epochTime = (epochStart:dt:epochEnd)';

        baseMask = epochTime >= baselineWindowSec(1) & epochTime <= baselineWindowSec(2);
        blockMask = epochTime >= 0 & epochTime <= condList(ci).blockDurationSecUsed;

        blockCount = condList(ci).blockCount;

        results.task.name = condName;
        results.task.blockOnsets = condList(ci).blockOnsets;
        results.task.blockDurationSec = condList(ci).blockDurationSecUsed;
        results.task.blockDurations = condList(ci).blockDurations;
        results.task.blockStartIdx = condList(ci).blockStartIdx;
        results.task.blockGapSec = blockGapSec;
        results.task.blockCount = blockCount;
        if isfield(condList(ci), 'blockLabels')
            results.task.blockLabels = condList(ci).blockLabels;
        end

        for ri = 1:nRoi
            hboMat = nan(numel(epochTime), blockCount);
            hbrMat = nan(numel(epochTime), blockCount);
            hboMeans = nan(blockCount,1);
            hbrMeans = nan(blockCount,1);

            for bi = 1:blockCount
                tQuery = condList(ci).blockOnsets(bi) + epochTime;
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

            results.rois(ri).name = roiDefs(ri).name;
            results.rois(ri).epochTime = epochTime;
            results.rois(ri).hboAvg = mean(hboMat, 2, 'omitnan');
            results.rois(ri).hbrAvg = mean(hbrMat, 2, 'omitnan');
            results.rois(ri).hboSem = hboSem;
            results.rois(ri).hbrSem = hbrSem;
            results.rois(ri).hboCi = hboCi;
            results.rois(ri).hbrCi = hbrCi;
            results.rois(ri).hboBlockMeans = hboMeans;
            results.rois(ri).hbrBlockMeans = hbrMeans;
            results.rois(ri).nHboChannels = roiNhbo(ri);
            results.rois(ri).nHbrChannels = roiNhbr(ri);

            for bi = 1:blockCount
                blockLabel = '';
                if isfield(condList(ci), 'blockLabels') && numel(condList(ci).blockLabels) >= bi
                    blockLabel = condList(ci).blockLabels{bi};
                end
                rowCell(end+1,:) = {subName, condName, roiDefs(ri).name, 'HbO', bi, blockLabel, hboMeans(bi), roiNhbo(ri)}; %#ok<AGROW>
                rowCell(end+1,:) = {subName, condName, roiDefs(ri).name, 'HbR', bi, blockLabel, hbrMeans(bi), roiNhbr(ri)}; %#ok<AGROW>
            end
        end
    end

    %% ===================== [8] SAVE OUTPUTS =====================
    outMat = fullfile(outDir, sprintf('%s_task-cate_blockavg.mat', subName));
    save(outMat, 'results', 'stim', 'roiDefs', '-v7.3');

    outCsv = fullfile(outDir, sprintf('%s_task-cate_blockavg_table.csv', subName));
    T = cell2table(rowCell, 'VariableNames', ...
        {'subjId','condition','roi','chromophore','blockIndex','blockLabel','blockMean','nChannels'});
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
    figPath0 = fullfile(figDir, sprintf('%s_task-cate_step-blockavg_qc_blockboundaries.png', subName));
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
    figPath1 = fullfile(figDir, sprintf('%s_task-cate_step-blockavg_qc_channels.png', subName));
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
    figPath2 = fullfile(figDir, sprintf('%s_task-cate_step-blockavg_signal.png', subName));
    saveas(fig2, figPath2);
    close(fig2);

    % Result plot: block average time course
    fig3 = figure('Visible','off');
    for ri = 1:nRoi
        subplot(nRoi,1,ri);
        x = results.rois(ri).epochTime;
        hboAvg = results.rois(ri).hboAvg;
        hbrAvg = results.rois(ri).hbrAvg;
        hboCi = results.rois(ri).hboCi;
        hbrCi = results.rois(ri).hbrCi;
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
        ylabel(sprintf('%s | %s', results.task.name, results.rois(ri).name), 'Interpreter','none');
        if ri == 1
            title(sprintf('%s block average (baseline-corrected)', subName), 'Interpreter','none');
        end
        if ri == nRoi
            xlabel('Time (s)');
        end
    end
    figPath3 = fullfile(figDir, sprintf('%s_task-cate_step-blockavg_result.png', subName));
    saveas(fig3, figPath3);
    close(fig3);

    % Trend plot: block index vs block mean (drift/fatigue check)
    fig4 = figure('Visible','off');
    for ri = 1:nRoi
        subplot(nRoi,1,ri);
        x = 1:numel(results.rois(ri).hboBlockMeans);
        plot(x, results.rois(ri).hboBlockMeans, '-o', 'Color', colorHbO, 'LineWidth', 1); hold on;
        plot(x, results.rois(ri).hbrBlockMeans, '-o', 'Color', colorHbR, 'LineWidth', 1);
        ylabel(sprintf('%s | %s', results.task.name, results.rois(ri).name), 'Interpreter','none');
        if ri == 1
            title(sprintf('%s block mean trend', subName), 'Interpreter','none');
            legend({'HbO','HbR'}, 'Location','best');
        end
        if ri == nRoi
            xlabel('Block index');
        end
        grid on;
    end
    figPath4 = fullfile(figDir, sprintf('%s_task-cate_step-blockavg_trend.png', subName));
    saveas(fig4, figPath4);
    close(fig4);

    fprintf('[DONE] %s | blocks=%d | out=%s\n', subName, ...
        condList(1).blockCount, outMat);

    if processOnlyFirst
        break;
    end
end
