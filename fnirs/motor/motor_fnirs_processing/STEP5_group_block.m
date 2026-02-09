%% STEP5_group_block.m
% Group-level analysis for block-average outputs (STEP3_BlockAvarage.m)
%
% Assumptions & Notes:
% - Inputs are per-subject blockavg tables: sub-xxx/blockavg/*_blockavg_table.csv
% - Tables include "condition" (left/right) and ROI (RightCentral/LeftCentral).
% - Group stats are one-sample tests vs 0 on per-subject mean block responses.
% - Optional group time-course averages are computed from *_blockavg.mat.

clear; clc;

%% ===================== [1] CONFIG =====================
% Root derived from this script location
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
rootDir = fileparts(thisDir);

% Subjects / tasks
subjList = {};               % empty = auto-discover sub-*
taskList = {'motor'};         % used for input/output naming

% Input naming
blockavgDirName = 'blockavg';
tableSuffixFmt = '_task-%s_blockavg_table.csv';
matSuffixFmt = '_task-%s_blockavg.mat';

% Output directories
outDir = fullfile(thisDir,  'group');
figDir = fullfile(thisDir,  'group');

% Stats options
doFdr = true;
fdrAlpha = 0.05;

% Plot options
condOrder = {'left','right'};
roiOrder = {'RightCentral','LeftCentral'};
chromOrder = {'HbO','HbR'};
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue

%% ===================== [2] DISCOVER SUBJECTS =====================
if isempty(subjList)
    subDirs = dir(fullfile(rootDir, 'sub-*'));
    if isempty(subDirs)
        error('No sub-* folders found under: %s', rootDir);
    end
    subjList = {subDirs.name};
end

if ~exist(outDir, 'dir'); mkdir(outDir); end
if ~exist(figDir, 'dir'); mkdir(figDir); end
fprintf('[INFO] Assumptions: condition=left/right, ROI=RightCentral/LeftCentral\n');

%% ===================== [3] LOAD SUBJECT BLOCKAVG TABLES =====================
allTables = {};
for si = 1:numel(subjList)
    subName = subjList{si};
    subPath = fullfile(rootDir, subName);
    for ti = 1:numel(taskList)
        taskName = taskList{ti};
        inCsv = fullfile(subPath, blockavgDirName, sprintf([subName tableSuffixFmt], taskName));
        if exist(inCsv, 'file') ~= 2
            % fallback: pick any *_blockavg_table.csv in blockavg dir
            tmp = dir(fullfile(subPath, blockavgDirName, '*_blockavg_table.csv'));
            if isempty(tmp)
                fprintf('[SKIP] %s %s: blockavg table not found.\n', subName, taskName);
                continue;
            end
            inCsv = fullfile(tmp(1).folder, tmp(1).name);
            tok = regexp(tmp(1).name, '_task-([^_]+)_blockavg_table\.csv', 'tokens', 'once');
            if ~isempty(tok)
                taskName = tok{1};
            end
        end

        T = readtable(inCsv);
        reqVars = {'subjId','condition','roi','chromophore','blockIndex','blockMean','nChannels'};
        if ~all(ismember(reqVars, T.Properties.VariableNames))
            error('Missing required columns in %s', inCsv);
        end
        T.subjId = string(T.subjId);
        T.condition = string(T.condition);
        T.roi = string(T.roi);
        T.chromophore = string(T.chromophore);
        if ~ismember('task', T.Properties.VariableNames)
            T.task = repmat(string(taskName), height(T), 1);
        else
            T.task = string(T.task);
        end
        allTables{end+1} = T; %#ok<AGROW>
    end
end

if isempty(allTables)
    error('No blockavg tables loaded. Check STEP3_BlockAvarage outputs.');
end

allT = vertcat(allTables{:});
condFromData = cellstr(unique(allT.condition));
condList = condOrder(ismember(condOrder, condFromData));
extraConds = setdiff(condFromData, condOrder, 'stable');
condList = [condList, extraConds];
if isempty(condList)
    condList = condFromData;
end

%% ===================== [4] SUBJECT-LEVEL SUMMARY =====================
[g, gSub, gTask, gCond, gRoi, gChrom] = findgroups(allT.subjId, allT.task, allT.condition, allT.roi, allT.chromophore);

blockMean = allT.blockMean;
nChannels = allT.nChannels;

meanBlock = splitapply(@(x) mean(x, 'omitnan'), blockMean, g);
sdBlock = splitapply(@(x) std(x, 'omitnan'), blockMean, g);
nBlocks = splitapply(@(x) sum(isfinite(x)), blockMean, g);
seBlock = sdBlock ./ sqrt(max(nBlocks,1));
meanNChannels = splitapply(@(x) mean(x, 'omitnan'), nChannels, g);

subT = table(gSub, gTask, gCond, gRoi, gChrom, nBlocks, meanBlock, sdBlock, seBlock, meanNChannels, ...
    'VariableNames', {'subjId','task','condition','roi','chromophore','nBlocks','meanBlock','sdBlock','seBlock','meanNChannels'});

outSubCsv = fullfile(outDir, 'group_blockavg_roi_subjects.csv');
writetable(subT, outSubCsv);

%% ===================== [5] GROUP STATS =====================
tasks = unique(subT.task);
condListUse = condList;
roiList = unique(subT.roi);
chromList = unique(subT.chromophore);

rowCell = {};
for ti = 1:numel(tasks)
    taskName = tasks(ti);
    for ci2 = 1:numel(condListUse)
        condName = string(condListUse{ci2});
        for ri = 1:numel(roiList)
            roiName = roiList(ri);
            for ci = 1:numel(chromList)
                chromName = chromList(ci);
                mask = (subT.task == taskName) & (subT.condition == condName) & ...
                    (subT.roi == roiName) & (subT.chromophore == chromName);
                vals = subT.meanBlock(mask);
                vals = vals(isfinite(vals));
                nSub = numel(vals);
                if nSub > 0
                    meanVal = mean(vals, 'omitnan');
                    sdVal = std(vals, 'omitnan');
                    seVal = sdVal / sqrt(nSub);
                else
                    meanVal = NaN; sdVal = NaN; seVal = NaN;
                end
                if nSub >= 2 && seVal > 0
                    tVal = meanVal / seVal;
                    dof = nSub - 1;
                    if exist('tcdf', 'file') == 2
                        pVal = 2 * tcdf(-abs(tVal), dof);
                    else
                        pVal = NaN;
                    end
                else
                    tVal = NaN; pVal = NaN; dof = nSub - 1;
                end
                if nSub >= 2 && sdVal > 0
                    cohend = meanVal / sdVal;
                else
                    cohend = NaN;
                end
                rowCell(end+1,:) = {string(taskName), string(condName), string(roiName), string(chromName), ...
                    nSub, meanVal, sdVal, seVal, tVal, pVal, dof, cohend}; %#ok<AGROW>
            end
        end
    end
end

groupT = cell2table(rowCell, 'VariableNames', ...
    {'task','condition','roi','chromophore','nSub','meanBlock','sdBlock','seBlock','t','p','dof','cohenD'});
groupT.task = string(groupT.task);
groupT.condition = string(groupT.condition);
groupT.roi = string(groupT.roi);
groupT.chromophore = string(groupT.chromophore);

% FDR correction (Benjamini-Hochberg)
if doFdr
    p = groupT.p;
    mask = isfinite(p);
    pvec = p(mask);
    qvec = nan(size(pvec));
    if ~isempty(pvec)
        [ps, idx] = sort(pvec);
        m = numel(ps);
        for i = 1:m
            qvec(i) = ps(i) * m / i;
        end
        for i = m-1:-1:1
            qvec(i) = min(qvec(i), qvec(i+1));
        end
        q = nan(size(p));
        qtmp = nan(size(pvec));
        qtmp(idx) = qvec; % map sorted q back to original p order
        q(mask) = qtmp;
        groupT.q = q;
        groupT.qSig = groupT.q < fdrAlpha;
    else
        groupT.q = nan(size(p));
        groupT.qSig = false(size(p));
    end
else
    groupT.q = nan(height(groupT),1);
    groupT.qSig = false(height(groupT),1);
end

outGroupCsv = fullfile(outDir, 'group_blockavg_roi_stats.csv');
writetable(groupT, outGroupCsv);

fprintf('[DONE] Group tables saved: %s\n', outGroupCsv);

%% ===================== [6] GROUP TIME COURSE (OPTIONAL) =====================
groupTC = struct();
for ti = 1:numel(tasks)
    taskName = tasks(ti);
    groupTC(ti).task = char(taskName);
    groupTC(ti).condition = struct();
    condCount = 0;

    for ci = 1:numel(condListUse)
        condName = string(condListUse{ci});
        epochTimeBase = [];
        roiName = "";
        hboAll = [];
        hbrAll = [];

        for si = 1:numel(subjList)
            subName = subjList{si};
            subPath = fullfile(rootDir, subName);
            inMat = fullfile(subPath, blockavgDirName, sprintf([subName matSuffixFmt], taskName));
            if exist(inMat, 'file') ~= 2
                tmp = dir(fullfile(subPath, blockavgDirName, '*_blockavg.mat'));
                if isempty(tmp)
                    continue;
                end
                inMat = fullfile(tmp(1).folder, tmp(1).name);
            end

            S = load(inMat, 'results');
            if ~isfield(S, 'results') || ~isfield(S.results, 'conditions')
                continue;
            end
            res = S.results;
            if isempty(res.conditions)
                continue;
            end
            condNames = string({res.conditions.name});
            idx = find(condNames == condName, 1);
            if isempty(idx)
                continue;
            end
            cRes = res.conditions(idx);

            if isempty(epochTimeBase)
                epochTimeBase = cRes.epochTime(:);
                roiName = string(cRes.roiName);
            end

            hboVec = cRes.hboAvg(:);
            hbrVec = cRes.hbrAvg(:);

            if numel(cRes.epochTime) ~= numel(epochTimeBase) || ...
                    max(abs(cRes.epochTime(:) - epochTimeBase)) > 1e-6
                hboVec = interp1(cRes.epochTime(:), hboVec, epochTimeBase, 'linear', NaN);
                hbrVec = interp1(cRes.epochTime(:), hbrVec, epochTimeBase, 'linear', NaN);
            end

            hboAll = [hboAll, hboVec]; %#ok<AGROW>
            hbrAll = [hbrAll, hbrVec]; %#ok<AGROW>
        end

        if isempty(epochTimeBase)
            continue;
        end

        hboValid = any(isfinite(hboAll), 1);
        hbrValid = any(isfinite(hbrAll), 1);
        hboMat = hboAll(:, hboValid);
        hbrMat = hbrAll(:, hbrValid);

        nSubHbo = size(hboMat, 2);
        nSubHbr = size(hbrMat, 2);

        hboMean = mean(hboMat, 2, 'omitnan');
        hbrMean = mean(hbrMat, 2, 'omitnan');

        hboSd = std(hboMat, 0, 2, 'omitnan');
        hbrSd = std(hbrMat, 0, 2, 'omitnan');

        hboSe = hboSd / sqrt(max(nSubHbo,1));
        hbrSe = hbrSd / sqrt(max(nSubHbr,1));

        condCount = condCount + 1;
        groupTC(ti).condition(condCount).name = char(condName);
        groupTC(ti).condition(condCount).roiName = char(roiName);
        groupTC(ti).condition(condCount).epochTime = epochTimeBase;
        groupTC(ti).condition(condCount).hboMean = hboMean;
        groupTC(ti).condition(condCount).hbrMean = hbrMean;
        groupTC(ti).condition(condCount).hboSe = hboSe;
        groupTC(ti).condition(condCount).hbrSe = hbrSe;
        groupTC(ti).condition(condCount).nSubHbo = nSubHbo;
        groupTC(ti).condition(condCount).nSubHbr = nSubHbr;
    end

    if condCount == 0
        continue;
    end

    outMat = fullfile(outDir, sprintf('group_task-%s_blockavg_timecourse.mat', taskName));
    save(outMat, 'groupTC', '-v7.3');

    % Plot time course (mean +/- SE)
    fig1 = figure('Visible','off');
    nCond = condCount;
    for ci = 1:nCond
        subplot(nCond,1,ci);
        t = groupTC(ti).condition(ci).epochTime;

        hboMean = groupTC(ti).condition(ci).hboMean;
        hbrMean = groupTC(ti).condition(ci).hbrMean;
        hboSe = groupTC(ti).condition(ci).hboSe;
        hbrSe = groupTC(ti).condition(ci).hbrSe;

        if any(isfinite(hboSe))
            fill([t; flipud(t)], [hboMean - hboSe; flipud(hboMean + hboSe)], ...
                colorHbO, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;
        end
        if any(isfinite(hbrSe))
            fill([t; flipud(t)], [hbrMean - hbrSe; flipud(hbrMean + hbrSe)], ...
                colorHbR, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;
        end

        plot(t, hboMean, 'Color', colorHbO, 'LineWidth', 1.5); hold on;
        plot(t, hbrMean, 'Color', colorHbR, 'LineWidth', 1.5);
        xline(0, '--k');
        ylabel(sprintf('%s | %s', groupTC(ti).condition(ci).name, ...
            groupTC(ti).condition(ci).roiName), 'Interpreter','none');
        if ci == 1
            title(sprintf('%s group block average', taskName), 'Interpreter','none');
        end
        if ci == nCond
            xlabel('Time (s)');
        end
    end
    figPath = fullfile(figDir, sprintf('group_task-%s_step-blockavg_timecourse.png', taskName));
    saveas(fig1, figPath);
    close(fig1);
end

%% ===================== [7] GROUP BAR PLOT =====================
for ti = 1:numel(tasks)
    taskName = tasks(ti);
    for ci = 1:numel(condListUse)
        condName = string(condListUse{ci});
        mask = (groupT.task == taskName) & (groupT.condition == condName);
        gPlotT = groupT(mask,:);

        roiNames = roiOrder(ismember(roiOrder, cellstr(unique(gPlotT.roi))));
        chromNames = chromOrder(ismember(chromOrder, cellstr(unique(gPlotT.chromophore))));

        if isempty(roiNames) || isempty(chromNames)
            continue;
        end

        data = nan(numel(roiNames), numel(chromNames));
        se = nan(numel(roiNames), numel(chromNames));
        for ri = 1:numel(roiNames)
            for ci2 = 1:numel(chromNames)
                idx = (gPlotT.roi == roiNames{ri}) & (gPlotT.chromophore == chromNames{ci2});
                if any(idx)
                    data(ri,ci2) = gPlotT.meanBlock(idx);
                    se(ri,ci2) = gPlotT.seBlock(idx);
                end
            end
        end

        fig1 = figure('Visible','off');
        b = bar(categorical(roiNames), data);
        if numel(b) >= 1
            b(1).FaceColor = colorHbO;
        end
        if numel(b) >= 2
            b(2).FaceColor = colorHbR;
        end
        hold on;
        for bi = 1:numel(b)
            if isprop(b(bi), 'XEndPoints')
                x = b(bi).XEndPoints;
            else
                x = b(bi).XData + b(bi).XOffset;
            end
            errorbar(x, data(:,bi), se(:,bi), 'k', 'LineStyle', 'none', 'LineWidth', 1);
        end
        ylabel('Block mean');
        legend(chromNames, 'Location','best');
        title(sprintf('%s %s group ROI block mean', taskName, condName), 'Interpreter','none');
        figPath = fullfile(figDir, sprintf('group_task-%s_cond-%s_step-blockavg_roi_mean.png', ...
            taskName, condName));
        saveas(fig1, figPath);
        close(fig1);
    end
end
