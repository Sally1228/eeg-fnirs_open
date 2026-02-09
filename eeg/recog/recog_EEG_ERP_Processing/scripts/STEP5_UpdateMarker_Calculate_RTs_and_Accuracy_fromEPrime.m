%Script #5 (E-Prime trials + EEG event integration)
%Extract trial-level behavior from E-Prime logs (memory/recognition) and write into EEG events.
% 重编码规则
% 记忆阶段，marker=4，renspose=1 → 11，mem1
% 记忆阶段，marker=4，renspose=2 → 12，mem2
% 记忆阶段，marker=4，没有回答/或者反应时间为NaN等→1*（*的意思是数据缺失），mem_na
% 
% 再认阶段，marker=4，回答正确 → 211，recog_old_right
% 再认阶段，marker=4，回答错误 → 212，recog_old_wrong
% 再认阶段，marker=4，没有回答/或者反应时间为NaN等→ 21*，recog_old_na
% 
% 再认阶段，marker=8，回答正确 → 221，recog_new_right
% 再认阶段，marker=8，回答错误 → 222，recog_new_wrong
% 再认阶段，marker=8，没有回答/或者反应时间为NaN等→ 22*，recog_new_na

%% Config
close all; clearvars;

% Assumptions & Notes:
% - E-Prime logs are stored in inputs/Behavior_Measurements/eprime_data.
% - Subject folders are stored in the sibling data directory (eeg/recog/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

input_dir = fullfile(project_dir, 'inputs');
output_dir = fullfile(project_dir, 'outputs');
EPrime_Dir = fullfile(input_dir, 'eprime_data');
data_dir = fileparts(project_dir);

fprintf('[%s] Assumption: E-Prime logs in %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), EPrime_Dir);
fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP5 Update markers + E-Prime RTs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to process
%'sub-009' 没有eprime数据
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-008', 'sub-010'};  

task = 'recog';

%E-Prime procedures for memory/recognition trials
Mem_Procedure = 'membase';
Recog_Procedure = 'regtest';
Mem_Prac_Run_Keyword = 'prac';

%Optional: skip first N stimulus events in EEG (e.g., practice) EEG中没有训练数据？
Skip_Stim_Events = 0;

%Whether to replace EEG.event.type with response-coded markers
Update_EventType = true; % true -> 11/12/211/212/221/222, false -> keep 4/8

Trial_CSV = fullfile(output_dir, 'EPrime_Trials_All_Subjects_recog.csv');
Summary_CSV = fullfile(output_dir, 'EPrime_RTs_Trial_Counts_&_ACCs_recog.csv');

allTrialRows = {};
summaryRows = {};
allTrialRows(1, :) = {'SubID','Trial','Stage','Marker','Response','CorrectResp','ACC','RT','GroupCode','GroupLabel'};
summaryRows(1, :) = {'SubID','Mem1 Mean RT','Mem2 Mean RT','Mem1 Trial Count','Mem2 Trial Count', ...
    'Recog Same Correct RT','Recog Diff Correct RT','Recog Same Correct Count','Recog Diff Correct Count', ...
    'Recog Same Accuracy','Recog Diff Accuracy','Recog Total Accuracy'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));
    filePath = fullfile(EPrime_Dir, [SUB{i} '_' task '_eprime.txt']);
    Subject_Path = fullfile(data_dir, SUB{i});
    EEG_In_Name = [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr.set'];
    EEG_Out_Name = [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_marker.set'];
    EEG_In_File = fullfile(Subject_Path, EEG_In_Name);
    EEG_Out_File = fullfile(Subject_Path, EEG_Out_Name);
    if ~exist(filePath, 'file')
        warning('E-Prime file not found: %s', filePath);
        continue;
    end

    markerList = [];
    accList = [];
    rtList = [];
    respList = {};
    crespList = {};
    stageList = {};

    fid = fopen(filePath, 'r', 'n', 'UTF-16LE');
    if fid == -1
        warning('Cannot open E-Prime file: %s', filePath);
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
        warning('No test trials found in: %s', filePath);
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

    for t = 1:nTrials
        allTrialRows(end+1, :) = {SUB{i}, t, stageList{t}, markerList(t), respList{t}, ...
            crespList{t}, accList(t), rtList(t), groupCode{t}, groupLabel{t}}; %#ok<AGROW>
    end

    isMem = strcmp(stageList, 'mem');
    isRecog = strcmp(stageList, 'recog');

    noRespVec = isnan(respNumList) | isnan(rtList);
    mem1Idx = isMem & (respNumList == 1) & ~noRespVec;
    mem2Idx = isMem & (respNumList == 2) & ~noRespVec;

    mem1MeanRT = mean(rtList(mem1Idx), 'omitnan');
    mem2MeanRT = mean(rtList(mem2Idx), 'omitnan');
    mem1Count = sum(mem1Idx);
    mem2Count = sum(mem2Idx);

    sameIdx = isRecog & markerList == 4;
    diffIdx = isRecog & markerList == 8;
    correctIdx = accList == 1;

    sameCorrectIdx = sameIdx & correctIdx;
    diffCorrectIdx = diffIdx & correctIdx;

    meanSameRT = mean(rtList(sameCorrectIdx), 'omitnan');
    meanDiffRT = mean(rtList(diffCorrectIdx), 'omitnan');

    sameCorrectCount = sum(sameCorrectIdx);
    diffCorrectCount = sum(diffCorrectIdx);
    sameTotal = sum(sameIdx);
    diffTotal = sum(diffIdx);

    sameAcc = NaN;
    if sameTotal > 0
        sameAcc = sameCorrectCount / sameTotal * 100;
    end

    diffAcc = NaN;
    if diffTotal > 0
        diffAcc = diffCorrectCount / diffTotal * 100;
    end

    totalAcc = mean([sameAcc, diffAcc], 'omitnan');

    summaryRows(end+1, :) = {SUB{i}, mem1MeanRT, mem2MeanRT, mem1Count, mem2Count, ...
        meanSameRT, meanDiffRT, sameCorrectCount, diffCorrectCount, sameAcc, diffAcc, totalAcc}; %#ok<AGROW>

    if ~exist(EEG_In_File, 'file')
        warning('EEG file not found: %s', EEG_In_File);
        continue;
    end

    EEG = pop_loadset('filename', EEG_In_Name, 'filepath', Subject_Path);

    evtTypeStr = cell(1, numel(EEG.event));
    for k = 1:numel(EEG.event)
        tval = EEG.event(k).type;
        if isnumeric(tval)
            evtTypeStr{k} = num2str(tval);
        else
            evtTypeStr{k} = tval;
        end
    end
    evtTypeStr = cellfun(@strtrim, evtTypeStr, 'UniformOutput', false);

    %Remove leading test markers if present (pattern: 1 4 4 4 2)-代表测试练习
    markerMask = strcmp(evtTypeStr, '1') | strcmp(evtTypeStr, '2') | ...
        strcmp(evtTypeStr, '4') | strcmp(evtTypeStr, '8');
    markerIdx = find(markerMask);
    testPattern = {'1','4','4','4','2'};
    noTestPattern = {'1','4','4','4','4','4','4','4','4','4','4'};
    removedTestMarkers = false;
    if numel(markerIdx) >= 10 && ...
            isequal(evtTypeStr(markerIdx(1:numel(testPattern))), testPattern)
        EEG.event(markerIdx(1:10)) = [];
        removedTestMarkers = true;
        fprintf('Removed test markers for %s\n', SUB{i});
    elseif numel(markerIdx) >= numel(noTestPattern) && ...
            isequal(evtTypeStr(markerIdx(1:numel(noTestPattern))), noTestPattern)
        %No test markers to remove
    end

    if removedTestMarkers
        EEG = eeg_checkset(EEG, 'eventconsistency');
        evtTypeStr = cell(1, numel(EEG.event));
        for k = 1:numel(EEG.event)
            tval = EEG.event(k).type;
            if isnumeric(tval)
                evtTypeStr{k} = num2str(tval);
            else
                evtTypeStr{k} = tval;
            end
        end
        evtTypeStr = cellfun(@strtrim, evtTypeStr, 'UniformOutput', false);
    end

    isStim = strcmp(evtTypeStr, '4') | strcmp(evtTypeStr, '8');
    stimIdx = find(isStim);
    if Skip_Stim_Events > 0
        if Skip_Stim_Events >= numel(stimIdx)
            warning('Skip_Stim_Events too large for %s', SUB{i});
            continue;
        end
        stimIdx = stimIdx(Skip_Stim_Events+1:end);
    end

    if numel(stimIdx) ~= nTrials
        fprintf(2, 'Event count mismatch for %s: EEG=%d, E-Prime=%d\n', ...
            SUB{i}, numel(stimIdx), nTrials);
        error('Event count mismatch. Stopping for manual inspection.');
    end

    for t = 1:nTrials
        EEG.event(stimIdx(t)).eprime_marker = markerList(t);
        EEG.event(stimIdx(t)).eprime_resp = respList{t};
        EEG.event(stimIdx(t)).eprime_cresp = crespList{t};
        EEG.event(stimIdx(t)).eprime_acc = accList(t);
        EEG.event(stimIdx(t)).eprime_rt = rtList(t);
        EEG.event(stimIdx(t)).eprime_group = groupCode{t};
        EEG.event(stimIdx(t)).eprime_stage = stageList{t};

        if Update_EventType
            EEG.event(stimIdx(t)).type_raw = evtTypeStr{stimIdx(t)};
            if ~isempty(groupCode{t})
                EEG.event(stimIdx(t)).type = groupCode{t};
            end
        end
    end

    EEG = eeg_checkset(EEG, 'eventconsistency');
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_marker'],'savenew', [Subject_Path filesep SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_marker.set'],'overwrite','on', 'gui','off'); 

    EEG = pop_saveset(EEG, 'filename', EEG_Out_Name, 'filepath', Subject_Path);
end

writecell(allTrialRows, Trial_CSV);
writecell(summaryRows, Summary_CSV);
