%% UpdateMarker_fromEprime.m
%==========================================================================
% 本脚本用于依据 E-Prime 行为日志更新 fNIRS 预处理数据中的刺激标记（stim），
% 使其与 EEG 的重编码规则一致，便于后续耦联分析。
% 主要流程（编号）：
% 1) 输入/路径：定位 E-Prime 日志、fNIRS 预处理 MAT、被试列表与任务名。
% 2) 预处理步骤：读取 E-Prime 日志，提取记忆/再认试次，过滤练习试次。
% 3) 核心计算：按 EEG 规则重编码 (11/12/1*/211/212/21*/221/222/22*)，
%              并按时间顺序对齐 fNIRS 的 4/8 刺激事件。
% 4) 输出：将新 stim 写入 *_fnirs_preproc.mat（原地或带后缀保存）。
%==========================================================================

% Assumptions & Notes:
% - E-Prime 日志为 UTF-16LE 文本，文件名格式：sub-XXX_recog_eprime.txt。
% - fNIRS 预处理 MAT 中包含 stim（Homer3 StimClass 或 struct）。
% - fNIRS 中 4/8 事件顺序与 E-Prime 试次顺序一致（可选移除练习序列）。

clear; clc;

%% ===================== [1] CONFIG =====================
task = 'recog';

thisDir = fileparts(mfilename('fullpath'));
fnirsRootDir = fileparts(thisDir);               % .../fnirs/recog
projectDir = fileparts(fileparts(fnirsRootDir)); % repo root

eprimeDirCandidates = { ...
    fullfile(projectDir, 'eeg', 'recog', 'recog_EEG_ERP_Processing', 'inputs', 'eprime_data'), ...
    };

subjList = {};                  % empty -> auto-discover sub-*
includeBadQuality = false;      % include folders like "质量不好sub-010"
badQualityPrefix = '质量不好sub-';

preprocSuffix = '_task-recog_fnirs_preproc.mat';
outSuffix = '_task-recog_fnirs_preproc_marker.mat';
overwritePreproc = true;        % true: update in place; false: write new MAT with outSuffix

Mem_Procedure = 'membase';
Recog_Procedure = 'regtest';
Mem_Prac_Run_Keyword = 'prac';

Skip_Stim_Events = 0;           % optional: skip first N stim events (4/8)
RemoveTestPattern = true;       % remove leading practice pattern if detected
TestPattern = {'1','4','4','4','2'};
NoTestPattern = {'1','4','4','4','4','4','4','4','4','4','4'};
TestPatternRemoveCount = 10;    % follow EEG script behavior

CheckMarkerSequence = true;     % compare 4/8 sequence between E-Prime and fNIRS

stimDurationFallbackSec = 2;    % used if stim has only onset or onset+duration
stimValueFallback = 1;

OverwriteExistingGroupStim = true;
groupCodeOrder = {'11','12','1*','211','212','21*','221','222','22*'};

%% ===================== [2] RESOLVE PATHS & ASSUMPTIONS =====================
eprimeDir = '';
for ci = 1:numel(eprimeDirCandidates)
    if exist(eprimeDirCandidates{ci}, 'dir')
        eprimeDir = eprimeDirCandidates{ci};
        break;
    end
end
if isempty(eprimeDir)
    error('E-Prime directory not found. Checked: %s', strjoin(eprimeDirCandidates, ' | '));
end

fprintf('[%s] Assumption: E-Prime logs in %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), eprimeDir);
fprintf('[%s] Assumption: fNIRS preproc MAT contains stim with markers 4/8\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

if isempty(subjList)
    subDirs = dir(fullfile(fnirsRootDir, 'sub-*'));
    subjList = {subDirs.name};
    if includeBadQuality
        badDirs = dir(fullfile(fnirsRootDir, [badQualityPrefix '*']));
        subjList = [subjList, {badDirs.name}]; %#ok<AGROW>
    end
end
if isempty(subjList)
    error('No subject folders found under: %s', fnirsRootDir);
end

fprintf('[%s] Update markers from E-Prime (fNIRS)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ===================== [3] MAIN LOOP =====================
for i = 1:numel(subjList)
    subName = subjList{i};
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, i, numel(subjList));

    eprimeFile = fullfile(eprimeDir, [subName '_' task '_eprime.txt']);
    if ~exist(eprimeFile, 'file')
        warning('E-Prime file not found: %s', eprimeFile);
        continue;
    end

    subPath = fullfile(fnirsRootDir, subName);
    inMat = fullfile(subPath, [subName preprocSuffix]);
    if exist(inMat, 'file') ~= 2
        tmp = dir(fullfile(subPath, '*_preproc.mat'));
        if isempty(tmp)
            warning('No preproc MAT found for %s', subName);
            continue;
        end
        inMat = fullfile(tmp(1).folder, tmp(1).name);
    end

    %% ---- Parse E-Prime ----
    markerList = [];
    accList = [];
    rtList = [];
    respList = {};
    crespList = {};
    stageList = {};

    fid = fopen(eprimeFile, 'r', 'n', 'UTF-16LE');
    if fid == -1
        warning('Cannot open E-Prime file: %s', eprimeFile);
        continue;
    end

    inLogFrame = false;
    procedure = '';
    running = '';
    marker = NaN;
    stim1_acc = NaN;
    stim1_rt = NaN;
    stim1_resp = '';
    stim1_cresp = '';
    stim2_acc = NaN;
    stim2_rt = NaN;
    stim2_resp = '';
    stim2_cresp = '';

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~isempty(line) && line(1) == char(65279)
            line = line(2:end);
        end

        if strcmp(line, '*** LogFrame Start ***')
            inLogFrame = true;
            procedure = '';
            running = '';
            marker = NaN;
            stim1_acc = NaN;
            stim1_rt = NaN;
            stim1_resp = '';
            stim1_cresp = '';
            stim2_acc = NaN;
            stim2_rt = NaN;
            stim2_resp = '';
            stim2_cresp = '';
            continue;
        end

        if inLogFrame
            if startsWith(line, 'Procedure:')
                procedure = strtrim(extractAfter(line, 'Procedure:'));
            elseif startsWith(line, 'Running:')
                running = strtrim(extractAfter(line, 'Running:'));
            elseif startsWith(line, 'marker:')
                marker = str2double(extractAfter(line, 'marker:'));
            elseif startsWith(line, 'stim1.ACC:')
                stim1_acc = str2double(extractAfter(line, 'stim1.ACC:'));
            elseif startsWith(line, 'stim1.RT:')
                stim1_rt = str2double(extractAfter(line, 'stim1.RT:'));
            elseif startsWith(line, 'stim1.RESP:')
                stim1_resp = strtrim(extractAfter(line, 'stim1.RESP:'));
            elseif startsWith(line, 'stim1.CRESP:')
                stim1_cresp = strtrim(extractAfter(line, 'stim1.CRESP:'));
            elseif startsWith(line, 'stim2.ACC:')
                stim2_acc = str2double(extractAfter(line, 'stim2.ACC:'));
            elseif startsWith(line, 'stim2.RT:')
                stim2_rt = str2double(extractAfter(line, 'stim2.RT:'));
            elseif startsWith(line, 'stim2.RESP:')
                stim2_resp = strtrim(extractAfter(line, 'stim2.RESP:'));
            elseif startsWith(line, 'stim2.CRESP:')
                stim2_cresp = strtrim(extractAfter(line, 'stim2.CRESP:'));
            end
        end

        if strcmp(line, '*** LogFrame End ***')
            inLogFrame = false;

            isMemTrial = strcmp(procedure, Mem_Procedure) && ~contains(lower(running), Mem_Prac_Run_Keyword);
            isRecogTrial = strcmp(procedure, Recog_Procedure);

            if isMemTrial || isRecogTrial
                if isMemTrial
                    stage = 'mem';
                    acc = stim1_acc;
                    rt = stim1_rt;
                    resp = stim1_resp;
                    cresp = stim1_cresp;
                else
                    stage = 'recog';
                    acc = stim2_acc;
                    rt = stim2_rt;
                    resp = stim2_resp;
                    cresp = stim2_cresp;
                end

                if rt == 0
                    rt = NaN;
                end

                markerList(end+1, 1) = marker; %#ok<AGROW>
                accList(end+1, 1) = acc; %#ok<AGROW>
                rtList(end+1, 1) = rt; %#ok<AGROW>
                respList{end+1, 1} = resp; %#ok<AGROW>
                crespList{end+1, 1} = cresp; %#ok<AGROW>
                stageList{end+1, 1} = stage; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    nTrials = numel(accList);
    if nTrials == 0
        warning('No test trials found in: %s', eprimeFile);
        continue;
    end

    respNumList = NaN(nTrials, 1);
    for t = 1:nTrials
        if isnumeric(respList{t})
            respNumList(t) = respList{t};
        else
            respNumList(t) = str2double(respList{t});
        end
    end

    groupCode = repmat({''}, nTrials, 1);
    groupLabel = repmat({''}, nTrials, 1);
    for t = 1:nTrials
        noResp = isnan(respNumList(t)) || isnan(rtList(t));
        if strcmp(stageList{t}, 'mem')
            if markerList(t) == 4 && respNumList(t) == 1 && ~noResp
                groupCode{t} = '11';
                groupLabel{t} = 'mem1';
            elseif markerList(t) == 4 && respNumList(t) == 2 && ~noResp
                groupCode{t} = '12';
                groupLabel{t} = 'mem2';
            else
                groupCode{t} = '1*';
                groupLabel{t} = 'mem_na';
            end
        elseif strcmp(stageList{t}, 'recog')
            if markerList(t) == 4 && noResp
                groupCode{t} = '21*';
                groupLabel{t} = 'recog_old_na';
            elseif markerList(t) == 4 && accList(t) == 1
                groupCode{t} = '211';
                groupLabel{t} = 'recog_old_right';
            elseif markerList(t) == 4
                groupCode{t} = '212';
                groupLabel{t} = 'recog_old_wrong';
            elseif markerList(t) == 8 && noResp
                groupCode{t} = '22*';
                groupLabel{t} = 'recog_new_na';
            elseif markerList(t) == 8 && accList(t) == 1
                groupCode{t} = '221';
                groupLabel{t} = 'recog_new_right';
            elseif markerList(t) == 8
                groupCode{t} = '222';
                groupLabel{t} = 'recog_new_wrong';
            else
                groupCode{t} = '';
                groupLabel{t} = 'recog_na';
            end
        else
            groupCode{t} = '';
            groupLabel{t} = 'NA';
        end
    end

    %% ---- Load stim from preproc MAT ----
    stim = [];
    try
        matObj = matfile(inMat);
        vars = who(matObj);
        if ~ismember('stim', vars)
            error('Missing "stim" in %s', inMat);
        end
        stim = matObj.stim;
    catch ME
        warning('Failed to read stim with matfile: %s', ME.message);
        S = load(inMat, 'stim');
        if ~isfield(S, 'stim')
            error('Missing "stim" in %s', inMat);
        end
        stim = S.stim;
    end

    if iscell(stim)
        stim = [stim{:}];
    end
    if isempty(stim)
        warning('Empty stim for %s', subName);
        continue;
    end
    stim = stim(:);

    %% ---- Build event lists (markers 1/2/4/8) ----
    markerSet = {'1','2','4','8'};
    allMarkers = {};
    allMarkerNum = [];
    allOnsets = [];
    allDur = [];
    allAmp = [];
    allIsStim = [];

    stimMarkers = {};
    stimMarkerNum = [];
    stimOnsets = [];
    stimDur = [];
    stimAmp = [];

    for k = 1:numel(stim)
        if isprop(stim(k), 'name')
            stimName = stim(k).name;
        elseif isfield(stim(k), 'name')
            stimName = stim(k).name;
        else
            stimName = '';
        end

        if iscell(stimName) && numel(stimName) == 1
            stimName = stimName{1};
        end
        if isnumeric(stimName)
            stimNameStr = num2str(stimName);
        else
            stimNameStr = char(stimName);
        end
        stimNameStr = strtrim(stimNameStr);

        if ~any(strcmp(stimNameStr, markerSet))
            continue;
        end

        if isprop(stim(k), 'data')
            stimData = stim(k).data;
        elseif isfield(stim(k), 'data')
            stimData = stim(k).data;
        else
            stimData = [];
        end
        if isempty(stimData)
            continue;
        end

        stimData = double(stimData);
        if size(stimData, 2) == 1
            stimData = [stimData, ...
                repmat(stimDurationFallbackSec, size(stimData,1), 1), ...
                repmat(stimValueFallback, size(stimData,1), 1)];
        elseif size(stimData, 2) == 2
            stimData = [stimData, repmat(stimValueFallback, size(stimData,1), 1)];
        elseif size(stimData, 2) > 3
            stimData = stimData(:,1:3);
        end

        markerNum = str2double(stimNameStr);
        isStim = strcmp(stimNameStr, '4') || strcmp(stimNameStr, '8');

        for r = 1:size(stimData, 1)
            allMarkers{end+1, 1} = stimNameStr; %#ok<AGROW>
            allMarkerNum(end+1, 1) = markerNum; %#ok<AGROW>
            allOnsets(end+1, 1) = stimData(r, 1); %#ok<AGROW>
            allDur(end+1, 1) = stimData(r, 2); %#ok<AGROW>
            allAmp(end+1, 1) = stimData(r, 3); %#ok<AGROW>
            allIsStim(end+1, 1) = isStim; %#ok<AGROW>

            if isStim
                stimMarkers{end+1, 1} = stimNameStr; %#ok<AGROW>
                stimMarkerNum(end+1, 1) = markerNum; %#ok<AGROW>
                stimOnsets(end+1, 1) = stimData(r, 1); %#ok<AGROW>
                stimDur(end+1, 1) = stimData(r, 2); %#ok<AGROW>
                stimAmp(end+1, 1) = stimData(r, 3); %#ok<AGROW>
            end
        end
    end

    if isempty(allOnsets) || isempty(stimOnsets)
        warning('No stim markers (1/2/4/8) found for %s', subName);
        continue;
    end

    allMat = [allMarkerNum, allOnsets, allDur, allAmp];
    [~, allKeep] = unique(allMat, 'rows', 'stable');
    allMarkers = allMarkers(allKeep);
    allMarkerNum = allMarkerNum(allKeep);
    allOnsets = allOnsets(allKeep);
    allDur = allDur(allKeep);
    allAmp = allAmp(allKeep);
    allIsStim = allIsStim(allKeep);

    stimMat = [stimMarkerNum, stimOnsets, stimDur, stimAmp];
    [~, stimKeep] = unique(stimMat, 'rows', 'stable');
    stimMarkers = stimMarkers(stimKeep);
    stimMarkerNum = stimMarkerNum(stimKeep);
    stimOnsets = stimOnsets(stimKeep);
    stimDur = stimDur(stimKeep);
    stimAmp = stimAmp(stimKeep);

    [allOnsets, allOrder] = sort(allOnsets);
    allMarkers = allMarkers(allOrder);
    allMarkerNum = allMarkerNum(allOrder);
    allIsStim = allIsStim(allOrder);

    [stimOnsets, stimOrder] = sort(stimOnsets);
    stimMarkers = stimMarkers(stimOrder);
    stimMarkerNum = stimMarkerNum(stimOrder);
    stimDur = stimDur(stimOrder);
    stimAmp = stimAmp(stimOrder);

    removedStimCount = 0;
    if RemoveTestPattern && numel(allMarkers) >= numel(TestPattern)
        if isequal(allMarkers(1:numel(TestPattern)), TestPattern)
            removeCount = min(TestPatternRemoveCount, numel(allMarkers));
            removedStimCount = sum(allIsStim(1:removeCount));
            if removedStimCount > 0
                if removedStimCount >= numel(stimOnsets)
                    warning('Test pattern removal too large for %s', subName);
                    continue;
                end
                stimOnsets = stimOnsets(removedStimCount+1:end);
                stimMarkers = stimMarkers(removedStimCount+1:end);
                stimMarkerNum = stimMarkerNum(removedStimCount+1:end);
                stimDur = stimDur(removedStimCount+1:end);
                stimAmp = stimAmp(removedStimCount+1:end);
            end
            fprintf('[%s] Removed %d leading stim events for practice in %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), removedStimCount, subName);
        elseif numel(allMarkers) >= numel(NoTestPattern) && ...
                isequal(allMarkers(1:numel(NoTestPattern)), NoTestPattern)
            fprintf('[%s] No practice pattern removed for %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        end
    end

    if Skip_Stim_Events > 0
        if Skip_Stim_Events >= numel(stimOnsets)
            warning('Skip_Stim_Events too large for %s', subName);
            continue;
        end
        stimOnsets = stimOnsets(Skip_Stim_Events+1:end);
        stimMarkers = stimMarkers(Skip_Stim_Events+1:end);
        stimMarkerNum = stimMarkerNum(Skip_Stim_Events+1:end);
        stimDur = stimDur(Skip_Stim_Events+1:end);
        stimAmp = stimAmp(Skip_Stim_Events+1:end);
    end

    if numel(stimOnsets) ~= nTrials
        fprintf(2, 'Event count mismatch for %s: fNIRS=%d, E-Prime=%d\n', ...
            subName, numel(stimOnsets), nTrials);
        error('Event count mismatch. Stopping for manual inspection.');
    end

    if CheckMarkerSequence
        if any(isnan(markerList))
            warning('E-Prime marker list has NaN for %s', subName);
        else
            if ~isequal(markerList(:), stimMarkerNum(:))
                error('Marker sequence mismatch between E-Prime and fNIRS for %s', subName);
            end
        end
    end

    %% ---- Build new group stim entries ----
    stimUpdated = stim;
    if OverwriteExistingGroupStim
        keepMask = true(size(stimUpdated));
        for k = 1:numel(stimUpdated)
            if isprop(stimUpdated(k), 'name')
                nm = stimUpdated(k).name;
            elseif isfield(stimUpdated(k), 'name')
                nm = stimUpdated(k).name;
            else
                nm = '';
            end
            if iscell(nm) && numel(nm) == 1
                nm = nm{1};
            end
            if isnumeric(nm)
                nmStr = num2str(nm);
            else
                nmStr = char(nm);
            end
            nmStr = strtrim(nmStr);
            if any(strcmp(nmStr, groupCodeOrder))
                keepMask(k) = false;
            end
        end
        stimUpdated = stimUpdated(keepMask);
    end
    stimUpdated = stimUpdated(:);

    stimTemplate = stimUpdated(1);
    for g = 1:numel(groupCodeOrder)
        gcode = groupCodeOrder{g};
        idx = strcmp(groupCode, gcode);
        if ~any(idx)
            continue;
        end
        newData = [stimOnsets(idx), stimDur(idx), stimAmp(idx)];
        newStim = stimTemplate;
        if isprop(newStim, 'name')
            newStim.name = gcode;
        elseif isfield(newStim, 'name')
            newStim.name = gcode;
        end
        if isprop(newStim, 'data')
            newStim.data = newData;
        elseif isfield(newStim, 'data')
            newStim.data = newData;
        end
        if isprop(newStim, 'dataLabels')
            newStim.dataLabels = {'onset','duration','value'};
        elseif isfield(newStim, 'dataLabels')
            newStim.dataLabels = {'onset','duration','value'};
        end
        stimUpdated(end+1, 1) = newStim; %#ok<AGROW>
    end

    %% ---- Save updated stim ----
    if overwritePreproc
        outMat = inMat;
    else
        outMat = fullfile(subPath, [subName outSuffix]);
        if exist(outMat, 'file') ~= 2
            copyfile(inMat, outMat);
        end
    end

    try
        matObj = matfile(outMat, 'Writable', true);
        matObj.stim = stimUpdated;
    catch ME
        warning('matfile update failed for %s: %s', subName, ME.message);
        S = load(inMat);
        S.stim = stimUpdated;
        save(outMat, '-struct', 'S', '-v7.3');
    end

    fprintf('[%s] Updated stim saved: %s (trials=%d)\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), outMat, nTrials);
end
