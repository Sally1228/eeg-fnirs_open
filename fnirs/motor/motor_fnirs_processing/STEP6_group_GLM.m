%% STEP6_group_GLM.m
% Group-level analysis from subject ROI GLM tables (STEP4_GLM output)
%
% Assumptions & Notes:
% - Inputs are per-subject ROI summary CSVs: sub-xxx/glm/*_glm_roi_table.csv
% - Terms use left/right; ROIs are RightCentral/LeftCentral.
% - Group stats are one-sample tests vs 0 on ROI meanBeta values.
% - This is a second-level analysis (not a full mixed-effects model).

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
taskList = {'motor'};        % used for input/output naming
roiSuffixFmt = '_task-%s_glm_roi_table.csv';

% Output directories
outDir = fullfile(thisDir, 'group');
figDir = fullfile( thisDir, 'group');

% Analysis options
termsToUse = {'left','right'}; % empty = all terms
doFdr = true;
fdrAlpha = 0.05;
fdrIncludeChrom = true;     % true: family=(task,term) across ROI x chrom; false: family=(task,term,chrom)

% Plot options
plotTerms = ["left","right"];
roiOrder = {'RightCentral','LeftCentral'};
chromOrder = {'HbO','HbR'};
colorHbO = [1.00 0.00 0.00]; % red
colorHbR = [0.00 0.00 1.00]; % blue
plotShowQ = true;            % true: annotate q; false: annotate p
plotShowPAndQ = false;       % if true, show both (q first)

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
fprintf('[INFO] Assumptions: terms=left/right, ROI=RightCentral/LeftCentral\n');

%% ===================== [3] LOAD SUBJECT ROI TABLES =====================
allTables = {};
for si = 1:numel(subjList)
    subName = subjList{si};
    subPath = fullfile(rootDir, subName);
    for ti = 1:numel(taskList)
        taskName = taskList{ti};
        inCsv = fullfile(subPath, 'glm', sprintf([subName roiSuffixFmt], taskName));
        if exist(inCsv, 'file') ~= 2
            fprintf('[SKIP] %s %s: ROI table not found.\n', subName, taskName);
            continue;
        end
        T = readtable(inCsv);
        reqVars = {'subjId','task','roi','chromophore','term','meanBeta'};
        if ~all(ismember(reqVars, T.Properties.VariableNames))
            error('Missing required columns in %s', inCsv);
        end
        T.subjId = string(T.subjId);
        T.task = string(T.task);
        T.roi = string(T.roi);
        T.chromophore = string(T.chromophore);
        T.term = string(T.term);
        if ~isempty(termsToUse)
            termMask = ismember(T.term, string(termsToUse));
            T = T(termMask,:);
        end
        allTables{end+1} = T; %#ok<AGROW>
    end
end

if isempty(allTables)
    error('No ROI tables loaded. Check subject folders and STEP4_GLM outputs.');
end

allT = vertcat(allTables{:});

%% ===================== [4] GROUP STATS =====================
tasks = unique(allT.task);
roiList = unique(allT.roi);
chromList = unique(allT.chromophore);
termList = unique(allT.term);

rowCell = {};
for ti = 1:numel(tasks)
    taskName = tasks(ti);
    for ri = 1:numel(roiList)
        roiName = roiList(ri);
        for ci = 1:numel(chromList)
            chromName = chromList(ci);
            for ti2 = 1:numel(termList)
                termName = termList(ti2);
                mask = (allT.task == taskName) & (allT.roi == roiName) & ...
                       (allT.chromophore == chromName) & (allT.term == termName);
                vals = allT.meanBeta(mask);
                vals = vals(isfinite(vals));
                nSub = numel(vals);
                if nSub > 0
                    meanBeta = mean(vals, 'omitnan');
                    sdBeta = std(vals, 'omitnan');
                    seBeta = sdBeta / sqrt(nSub);
                else
                    meanBeta = NaN; sdBeta = NaN; seBeta = NaN;
                end
                if nSub >= 2 && seBeta > 0
                    tVal = meanBeta / seBeta;
                    dof = nSub - 1;
                    if exist('tcdf', 'file') == 2
                        pVal = 2 * tcdf(-abs(tVal), dof);
                    else
                        pVal = NaN;
                    end
                else
                    tVal = NaN; pVal = NaN; dof = nSub - 1;
                end
                if nSub >= 2 && sdBeta > 0
                    cohend = meanBeta / sdBeta;
                else
                    cohend = NaN;
                end
                rowCell(end+1,:) = {char(taskName), char(roiName), char(chromName), char(termName), ...
                    nSub, meanBeta, sdBeta, seBeta, tVal, pVal, dof, cohend}; %#ok<AGROW>
            end
        end
    end
end

groupT = cell2table(rowCell, 'VariableNames', ...
    {'task','roi','chromophore','term','nSub','meanBeta','sdBeta','seBeta','t','p','dof','cohenD'});
groupT.task = string(groupT.task);
groupT.roi = string(groupT.roi);
groupT.chromophore = string(groupT.chromophore);
groupT.term = string(groupT.term);

% FDR correction (Benjamini-Hochberg) within family
if doFdr
    groupT.q = nan(height(groupT),1);
    groupT.qSig = false(height(groupT),1);
    if fdrIncludeChrom
        famKey = groupT.task + "|" + groupT.term;
    else
        famKey = groupT.task + "|" + groupT.term + "|" + groupT.chromophore;
    end
    famList = unique(famKey);
    for fi = 1:numel(famList)
        idxFam = famKey == famList(fi);
        p = groupT.p(idxFam);
        mask = isfinite(p);
        if any(mask)
            pvec = p(mask);
            [ps, idx] = sort(pvec);
            m = numel(ps);
            qvec = nan(size(ps));
            for i = 1:m
                qvec(i) = ps(i) * m / i;
            end
            for i = m-1:-1:1
                qvec(i) = min(qvec(i), qvec(i+1));
            end
            q = nan(size(p));
            q(mask) = qvec(idx);
            groupT.q(idxFam) = q;
            groupT.qSig(idxFam) = q < fdrAlpha;
        end
    end
else
    groupT.q = nan(height(groupT),1);
    groupT.qSig = false(height(groupT),1);
end

outGroupCsv = fullfile(outDir, 'group_glm_roi_stats.csv');
writetable(groupT, outGroupCsv);

outSubCsv = fullfile(outDir, 'group_glm_roi_subjects.csv');
writetable(allT, outSubCsv);

fprintf('[DONE] Group tables saved: %s\n', outGroupCsv);

%% ===================== [5] FIGURES =====================
termCol = string(groupT.term);
taskCol = string(groupT.task);
plotTermsUse = plotTerms(ismember(plotTerms, unique(termCol)));
if isempty(plotTermsUse)
    plotTermsUse = unique(termCol);
end

for pti = 1:numel(plotTermsUse)
    plotTerm = string(plotTermsUse(pti));
    if ~any(termCol == plotTerm)
        continue;
    end
    for ti = 1:numel(tasks)
        taskName = tasks(ti);
        mask = (taskCol == taskName) & (termCol == plotTerm);
        subT = groupT(mask,:);

        subRoi = string(subT.roi);
        subChrom = string(subT.chromophore);
        roiNames = string(roiOrder);
        roiNames = roiNames(ismember(roiNames, unique(subRoi)));
        chromNames = string(chromOrder);
        chromNames = chromNames(ismember(chromNames, unique(subChrom)));

        if isempty(roiNames) || isempty(chromNames)
            continue;
        end

        data = nan(numel(roiNames), numel(chromNames));
        se = nan(numel(roiNames), numel(chromNames));
        pvals = nan(numel(roiNames), numel(chromNames));
        qvals = nan(numel(roiNames), numel(chromNames));
        for ri = 1:numel(roiNames)
            for ci = 1:numel(chromNames)
                idx = (subRoi == roiNames(ri)) & (subChrom == chromNames(ci));
                if any(idx)
                    rowIdx = find(idx, 1, 'first');
                    data(ri,ci) = subT.meanBeta(rowIdx);
                    se(ri,ci) = subT.seBeta(rowIdx);
                    pvals(ri,ci) = subT.p(rowIdx);
                    qvals(ri,ci) = subT.q(rowIdx);
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
        for ci = 1:numel(b)
            if isprop(b(ci), 'XEndPoints')
                x = b(ci).XEndPoints;
            else
                x = b(ci).XData + b(ci).XOffset;
            end
            errorbar(x, data(:,ci), se(:,ci), 'k', 'LineStyle', 'none', 'LineWidth', 1);
            for ri = 1:numel(roiNames)
                % subject-level points with jitter
                subjMask = (allT.task == taskName) & (allT.term == plotTerm) & ...
                    (allT.roi == roiNames(ri)) & (allT.chromophore == chromNames(ci));
                subjVals = allT.meanBeta(subjMask);
                subjVals = subjVals(isfinite(subjVals));
                if ~isempty(subjVals)
                    if exist('swarmchart', 'file') == 2
                        swarmchart(repmat(x(ri), numel(subjVals), 1), subjVals, 18, ...
                            'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', 'none', ...
                            'XJitter', 'density', 'XJitterWidth', 0.15);
                    else
                        jitter = (rand(numel(subjVals),1) - 0.5) * 0.15;
                        scatter(x(ri) + jitter, subjVals, 18, ...
                            'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', 'none');
                    end
                end
                if isfinite(pvals(ri,ci)) && isfinite(data(ri,ci)) && isfinite(se(ri,ci))
                    y = data(ri,ci) + se(ri,ci);
                    if plotShowPAndQ && isfinite(qvals(ri,ci))
                        label = sprintf('q=%.3f, p=%.3f', qvals(ri,ci), pvals(ri,ci));
                    elseif plotShowQ && isfinite(qvals(ri,ci))
                        label = sprintf('q=%.3f', qvals(ri,ci));
                    else
                        label = sprintf('p=%.3f', pvals(ri,ci));
                    end
                    text(x(ri), y, label, ...
                        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                        'FontSize', 8);
                end
            end
        end
        ylabel(sprintf('Beta (%s)', plotTerm));
        legend(chromNames, 'Location','best');
        title(sprintf('%s group ROI term: %s', taskName, plotTerm), 'Interpreter','none');
        figPath = fullfile(figDir, sprintf('group_task-%s_term-%s_step-glm_roi.png', taskName, plotTerm));
        saveas(fig1, figPath);
        close(fig1);
    end
end
