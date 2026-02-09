% EEG-fNIRS time-lag correlation and coherence analysis for auditory oddball task.
% Assumptions & Notes:
% - EEG bins: common=1, oddball=2 (update if your ERPLAB bin list differs).
% - fNIRS stim(1)=common, stim(2)=oddball in the preproc MAT.
% - EEG and fNIRS are aligned by event timing when possible; fallback is order-based.
% - fNIRS trial validity uses tIncAuto (global); channel-wise tIncAutoCh is not used.
% - EEG band power is computed from the continuous dataset and interpolated to fNIRS time.
% - Oddball-common coupling is computed as (oddball - common) on xcorr/coherence derived
%   from matched-trial masks (not by subtracting raw time series).


% 整体目标：以“连续EEG的频带功率包络”与“fNIRS HbO/HbR时序”做时滞相关与相干分析，同时显式处理EEG与fNIRS试次剔除不一致的问题。
% 核心流程
% - 配置区：定义EEG/ fNIRS数据路径模式、频段、ROI、对齐参数、相关与相干参数、输出路径等，确保脚本端到端可运行（配置在最上面）。
% - EEG事件与剔除：先加载“EEG分段数据”以读取 EEG.epoch.eventbini 和 EEG.reject，得到每个事件的类型（common/oddball）、时间戳、是否被剔除。
% - fNIRS事件与剔除：加载 dc/stim/tIncAuto，解析 stim(1) 与 stim(2) 的事件表，并用 tIncAuto 判断每个事件是否有效。
% 对齐策略：
% - 若EEG和fNIRS事件数量与类型一致，直接按序匹配。
% - 否则做“基于时间戳的匹配”，估计最优时间偏移 offsetSec，再按 alignToleranceSec 匹配事件。
% - 记录匹配结果与“双方都有效”的试次数量。
% EEG连续功率：
% - 载入“连续EEG”数据，选择通道（优先标签，如 Pz）。
% - 对每个频段带通滤波，计算Hilbert功率包络，并可选平滑。
% - 将EEG功率序列按 offsetSec 对齐后插值到fNIRS时间轴。
% fNIRS通道与ROI：
% - 读取 measurementList，区分HbO/HbR，依据 mlActAuto 剔除无效通道。
% - 可选ROI筛选（useRoiFilter）。
% 分析掩膜：
% - 默认 maskMode='matched_trials'：只在“EEG与fNIRS均有效且对齐的试次窗口”内做分析，窗口由 analysisWindowSec 控制。
% 跨相关与相干：
% - 对每个ROI × HbO/HbR × EEG频段，计算时滞相关峰值（xcorr）与相干谱（mscohere），输出峰值与频段统计。
% - Oddball-Common 使用 oddball 与 common 的 xcorr/coherence 差值来做对比。
% 输出：
% - STEP3_summary.csv：每名被试、Analysis、ROI、Hb类型、频段的峰值相关、滞后、相干统计。
% - STEP3_alignment_summary.csv：每名被试的EEG/fNIRS事件数、有效试次数、匹配数与估计偏移。
% 关键输出文件
% - STEP3_summary.csv
% - STEP3_alignment_summary.csv
% 关键参数建议关注
% - alignToleranceSec 和 minTimeMatchEvents：决定事件时间匹配的严格性。
% - maskMode 和 analysisWindowSec：决定只看试次窗口还是全时序。
% - maxLagSec、cohWindowSec、cohFreqRange：决定相关与相干分析的尺度。
% - useRoiFilter：是否按ROI分析。



%% Config
close all; clearvars; clc;

taskName = 'audi';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset patterns (continuous and epoched)
eegContPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins.set', ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist.set', ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr.set' ...
    };
eegEpochPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set', ...
    };

% EEG bin mapping
eegCommonBin = 1;
eegOddballBin = 2;

% EEG channels (labels preferred)
eegChanLabels = {'Pz'};
eegChanFallbackIdx = 53;
eegCombineMode = 'mean'; % 'mean' | 'median'

% EEG frequency bands
bandNames = {'Delta', 'Theta', 'Alpha', 'Beta'};
bandDefs = [1 4; 4 8; 8 13; 13 30];

% EEG power extraction
powerMethod = 'hilbert'; % 'hilbert'
smoothWindowSec = 0; % 0 = no smoothing

% fNIRS stim settings
stimIndexCommon = 1;
stimIndexOddball = 2;
stimDurationSec = 0.5;
fnirsTrialKeepFrac = 0.8;
fnirsTrialPadSec = 0.0;

% Alignment settings
alignToleranceSec = 0.3;
minTimeMatchEvents = 5;

% Analysis masks
maskMode = 'matched_trials'; % 'matched_trials' | 'all'
analysisWindowSec = [-2 20];

% Analysis conditions
analysisNames = {'All','Common','OddballCommon'};
analysisTags = {'all','common','oddballcommon'};
analysisTitles = {'All','Common','Oddball-Common'};

% Cross-correlation settings
maxLagSec = 10;
xcorrPeakMode = 'abs'; % 'abs' | 'pos'
normalizeMode = 'zscore'; % 'zscore' | 'demean' | 'none'

% Coherence settings
cohWindowSec = 20;
cohOverlapFrac = 0.5;
cohNfft = 256;
cohFreqRange = [0.01 0.2];

% fNIRS ROI settings (optional)
useRoiFilter = true;
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

% Output
outputDirName = 'STEP3_time-lag_correlation_and_coherence';
saveFullSpectra = false;

% Figure settings
saveFigures = true;
plotPerPairFigures = false;
plotSummaryFigures = true;
maxFiguresPerSubject = 12;
figureFormats = {'png'};
figureDpi = 150;
figureVisible = 'off'; % 'on' | 'off'
closeFigures = true;
summaryRoiName = 'All';
summaryChromList = {'HbO','HbR'};

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', outputDirName);
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end
figRoot = fullfile(outRoot, 'figures');
if saveFigures && exist(figRoot, 'dir') ~= 7
    mkdir(figRoot);
end

fprintf('[%s] Assumptions: EEG bins common=%d oddball=%d; stim(1)=common stim(2)=oddball; tIncAuto trial QC; timing alignment by events; oddball-common uses oddball minus common xcorr/coh.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), eegCommonBin, eegOddballBin);

if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %#ok<ASGLU>

if strcmpi(powerMethod, 'hilbert') && exist('hilbert', 'file') ~= 2
    error('hilbert not found (Signal Processing Toolbox required).');
end
if exist('mscohere', 'file') ~= 2
    error('mscohere not found (Signal Processing Toolbox required).');
end

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

%% Output tables
summaryHeader = {'SubID','Analysis','ROI','Chrom','Band','BandHz','SamplesUsed','PeakR','PeakLagSec', ...
    'CohMean','CohPeak','CohPeakHz','OffsetSec','MatchedEvents','MatchedValidEvents'};
summaryRows = cell(0, numel(summaryHeader));

alignHeader = {'SubID','EEGCommon','EEGCommonValid','EEGOdd','EEGOddValid', ...
    'FNIRSCommon','FNIRSCommonValid','FNIRSOdd','FNIRSOddValid', ...
    'MatchedEvents','MatchedValidEvents','OffsetSec'};
alignRows = cell(0, numel(alignHeader));

if numel(analysisNames) ~= numel(analysisTags) || numel(analysisNames) ~= numel(analysisTitles)
    error('analysisNames, analysisTags, analysisTitles must have the same length.');
end

%% Main loop
for si = 1:numel(subjList)
    subName = subjList{si};
    fprintf('\n[%s] Subject %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, si, numel(subjList));

    eegSubDir = fullfile(eegRoot, subName);
    fnirsSubDir = fullfile(fnirsRoot, subName);
    if exist(eegSubDir, 'dir') ~= 7 || exist(fnirsSubDir, 'dir') ~= 7
        fprintf('[%s] WARNING: Missing subject folder; skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end

    %% Locate EEG continuous dataset
    eegContPath = '';
    eegContFile = '';
    for pi = 1:numel(eegContPatterns)
        candidate = sprintf(eegContPatterns{pi}, subName, taskName);
        candidatePath = fullfile(eegSubDir, candidate);
        if exist(candidatePath, 'file') == 2
            eegContPath = candidatePath;
            eegContFile = candidate;
            break;
        end
    end
    if isempty(eegContPath)
        fprintf('[%s] WARNING: EEG continuous .set not found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    %% Locate EEG epoched dataset
    eegEpochPath = '';
    eegEpochFile = '';
    for pi = 1:numel(eegEpochPatterns)
        candidate = sprintf(eegEpochPatterns{pi}, subName, taskName);
        candidatePath = fullfile(eegSubDir, candidate);
        if exist(candidatePath, 'file') == 2
            eegEpochPath = candidatePath;
            eegEpochFile = candidate;
            break;
        end
    end
    if isempty(eegEpochPath)
        fprintf('[%s] WARNING: EEG epoched .set not found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    %% Locate fNIRS preproc MAT
    fnirsMat = fullfile(fnirsSubDir, sprintf('%s_task-%s_fnirs_preproc.mat', subName, taskName));
    if exist(fnirsMat, 'file') ~= 2
        fprintf('[%s] WARNING: fNIRS preproc MAT not found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    subOutDir = fullfile(outRoot, subName);
    if exist(subOutDir, 'dir') ~= 7
        mkdir(subOutDir);
    end
    figSubDir = fullfile(figRoot, subName);
    if saveFigures && plotPerPairFigures && exist(figSubDir, 'dir') ~= 7
        mkdir(figSubDir);
    end
    figCount = 0;

    %% Load EEG epoched dataset to get events and rejection
    EEG_ep = pop_loadset('filename', eegEpochFile, 'filepath', eegSubDir);
    if EEG_ep.trials < 1 || isempty(EEG_ep.epoch)
        fprintf('[%s] WARNING: EEG epoched data invalid for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    reject = false(1, EEG_ep.trials);
    if isfield(EEG_ep, 'reject') && ~isempty(EEG_ep.reject)
        rejectFields = fieldnames(EEG_ep.reject);
        for rf = 1:numel(rejectFields)
            val = EEG_ep.reject.(rejectFields{rf});
            if isempty(val) || ~(islogical(val) || isnumeric(val))
                continue;
            end
            if isvector(val) && numel(val) == EEG_ep.trials
                reject = reject | logical(val(:))';
            elseif ismatrix(val) && size(val, 2) == EEG_ep.trials
                reject = reject | any(val, 1);
            elseif ismatrix(val) && size(val, 1) == EEG_ep.trials
                reject = reject | any(val, 2)';
            end
        end
    end
    eegEventValid = ~reject(:);

    eegEventTimeSec = nan(EEG_ep.trials, 1);
    eegEventType = zeros(EEG_ep.trials, 1); % 1=common, 2=oddball

    for e = 1:EEG_ep.trials
        if ~isfield(EEG_ep.epoch(e), 'eventbini')
            continue;
        end
        eventBini = EEG_ep.epoch(e).eventbini;
        eventUrevent = [];
        if isfield(EEG_ep.epoch(e), 'eventurevent')
            eventUrevent = EEG_ep.epoch(e).eventurevent;
        end
        eventEvent = [];
        if isfield(EEG_ep.epoch(e), 'event')
            eventEvent = EEG_ep.epoch(e).event;
        end

        biniValsAll = [];
        if iscell(eventBini)
            for ev = 1:numel(eventBini)
                tmp = eventBini{ev};
                if iscell(tmp)
                    tmp = [tmp{:}];
                end
                if ischar(tmp)
                    tmp = str2double(tmp);
                end
                if isempty(tmp)
                    continue;
                end
                biniValsAll = [biniValsAll, double(tmp(:)')]; %#ok<AGROW>
            end
        else
            biniValsAll = double(eventBini(:)');
        end

        if any(biniValsAll == eegOddballBin)
            eegEventType(e) = 2;
            targetBin = eegOddballBin;
        elseif any(biniValsAll == eegCommonBin)
            eegEventType(e) = 1;
            targetBin = eegCommonBin;
        else
            continue;
        end

        lat = NaN;
        nEv = numel(eventBini);
        for ev = 1:nEv
            if iscell(eventBini)
                b = eventBini{ev};
            else
                b = eventBini(ev);
            end
            if iscell(b)
                b = [b{:}];
            end
            if ischar(b)
                b = str2double(b);
            end
            if isempty(b) || ~any(double(b) == targetBin)
                continue;
            end

            if ~isempty(eventUrevent) && isfield(EEG_ep, 'urevent') && ~isempty(EEG_ep.urevent)
                if iscell(eventUrevent)
                    u = eventUrevent{ev};
                else
                    u = eventUrevent(ev);
                end
                if iscell(u)
                    u = [u{:}];
                end
                if ~isempty(u) && u(1) >= 1 && u(1) <= numel(EEG_ep.urevent)
                    lat = EEG_ep.urevent(u(1)).latency;
                end
            end

            if ~isfinite(lat) && ~isempty(eventEvent) && isfield(EEG_ep, 'event') && ~isempty(EEG_ep.event)
                if iscell(eventEvent)
                    evIdx = eventEvent{ev};
                else
                    evIdx = eventEvent(ev);
                end
                if iscell(evIdx)
                    evIdx = [evIdx{:}];
                end
                if ~isempty(evIdx) && evIdx(1) >= 1 && evIdx(1) <= numel(EEG_ep.event)
                    lat = EEG_ep.event(evIdx(1)).latency;
                end
            end

            if isfinite(lat)
                break;
            end
        end

        if isfinite(lat)
            eegEventTimeSec(e) = (double(lat) - 1) / EEG_ep.srate;
        end
    end

    eegMask = eegEventType > 0;
    eegEventTimeSec = eegEventTimeSec(eegMask);
    eegEventType = eegEventType(eegMask);
    eegEventValid = eegEventValid(eegMask);

    eegCommonCount = sum(eegEventType == 1);
    eegOddCount = sum(eegEventType == 2);
    eegCommonValid = sum(eegEventType == 1 & eegEventValid);
    eegOddValid = sum(eegEventType == 2 & eegEventValid);

    %% Load fNIRS preproc
    S = load(fnirsMat, 'dc', 'stim', 'mlActAuto', 'tIncAuto', 'tIncAutoCh');
    if ~isfield(S, 'dc') || ~isfield(S, 'stim')
        fprintf('[%s] WARNING: fNIRS MAT missing dc/stim for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    dc = S.dc;
    stim = S.stim;
    if iscell(stim)
        stim = [stim{:}];
    end
    tFnirs = double(dc.time(:));
    yFnirs = double(dc.dataTimeSeries);
    if size(yFnirs,1) ~= numel(tFnirs)
        fprintf('[%s] WARNING: fNIRS time/data length mismatch for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    nStim = numel(stim);
    if nStim < max(stimIndexCommon, stimIndexOddball)
        fprintf('[%s] WARNING: stim has only %d entries for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nStim, subName);
        continue;
    end

    stimCommon = stim(stimIndexCommon).data;
    stimOddball = stim(stimIndexOddball).data;
    if isempty(stimCommon) || isempty(stimOddball)
        fprintf('[%s] WARNING: stim common/oddball empty for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    if size(stimCommon,2) == 1
        stimCommon = [stimCommon, repmat(stimDurationSec, size(stimCommon,1), 1), ones(size(stimCommon,1),1)];
    elseif size(stimCommon,2) == 2
        stimCommon = [stimCommon(:,1), stimCommon(:,2), ones(size(stimCommon,1),1)];
    else
        stimCommon = stimCommon(:,1:3);
    end
    if size(stimOddball,2) == 1
        stimOddball = [stimOddball, repmat(stimDurationSec, size(stimOddball,1), 1), ones(size(stimOddball,1),1)];
    elseif size(stimOddball,2) == 2
        stimOddball = [stimOddball(:,1), stimOddball(:,2), ones(size(stimOddball,1),1)];
    else
        stimOddball = stimOddball(:,1:3);
    end

    stimCommon = sortrows(stimCommon, 1);
    stimOddball = sortrows(stimOddball, 1);

    nCommonFnirs = size(stimCommon,1);
    nOddFnirs = size(stimOddball,1);

    fnirsEventType = [ones(nCommonFnirs,1); 2*ones(nOddFnirs,1)];
    fnirsEventTimeSec = [stimCommon(:,1); stimOddball(:,1)];
    fnirsEventDurSec = [stimCommon(:,2); stimOddball(:,2)];
    [fnirsEventTimeSec, sortIdx] = sort(fnirsEventTimeSec);
    fnirsEventType = fnirsEventType(sortIdx);
    fnirsEventDurSec = fnirsEventDurSec(sortIdx);

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
        if numel(tIncVec) ~= numel(tFnirs)
            if numel(tIncVec) > numel(tFnirs)
                tIncVec = tIncVec(1:numel(tFnirs));
            else
                tIncVec = [tIncVec; ones(numel(tFnirs) - numel(tIncVec), 1)];
            end
        end
    end

    fnirsEventValid = true(numel(fnirsEventTimeSec), 1);
    if ~isempty(tIncVec)
        for ei = 1:numel(fnirsEventTimeSec)
            onset = fnirsEventTimeSec(ei);
            dur = fnirsEventDurSec(ei);
            if ~isfinite(dur) || dur <= 0
                dur = stimDurationSec;
            end
            winStart = onset - fnirsTrialPadSec;
            winStop = onset + dur + fnirsTrialPadSec;
            idx = tFnirs >= winStart & tFnirs <= winStop;
            if any(idx)
                fnirsEventValid(ei) = mean(tIncVec(idx)) >= fnirsTrialKeepFrac;
            else
                fnirsEventValid(ei) = false;
            end
        end
    end

    fnirsCommonValid = sum(fnirsEventType == 1 & fnirsEventValid);
    fnirsOddValid = sum(fnirsEventType == 2 & fnirsEventValid);

    %% Align EEG and fNIRS events
    nEegEvents = numel(eegEventType);
    nFnirsEvents = numel(fnirsEventType);
    mapE = [];
    mapF = [];
    offsetSec = NaN;

    if nEegEvents == nFnirsEvents && all(eegEventType == fnirsEventType)
        mapE = (1:nEegEvents)';
        mapF = (1:nFnirsEvents)';
        if any(isfinite(eegEventTimeSec))
            dt = eegEventTimeSec - fnirsEventTimeSec;
            offsetSec = median(dt(isfinite(dt)));
        else
            offsetSec = 0;
        end
    else
        fprintf('[%s] WARNING: Event count/type mismatch (EEG=%d, fNIRS=%d). Using time-based matching.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nEegEvents, nFnirsEvents);

        nFinite = sum(isfinite(eegEventTimeSec));
        if nFinite >= minTimeMatchEvents
            offsets = [];
            nCandE = min(nEegEvents, 5);
            nCandF = min(nFnirsEvents, 5);
            for i = 1:nCandE
                if ~isfinite(eegEventTimeSec(i))
                    continue;
                end
                for j = 1:nCandF
                    if eegEventType(i) ~= fnirsEventType(j)
                        continue;
                    end
                    offsets(end+1) = eegEventTimeSec(i) - fnirsEventTimeSec(j); %#ok<AGROW>
                end
            end
            if isempty(offsets)
                offsets = eegEventTimeSec(isfinite(eegEventTimeSec(1:nCandE))) - fnirsEventTimeSec(1);
            end

            bestOffset = offsets(1);
            bestMatch = -inf;
            for oc = 1:numel(offsets)
                offset = offsets(oc);
                eIdx = 1;
                fIdx = 1;
                count = 0;
                while eIdx <= nEegEvents && fIdx <= nFnirsEvents
                    if ~isfinite(eegEventTimeSec(eIdx))
                        eIdx = eIdx + 1;
                        continue;
                    end
                    if eegEventType(eIdx) ~= fnirsEventType(fIdx)
                        dtType = (eegEventTimeSec(eIdx) - offset) - fnirsEventTimeSec(fIdx);
                        if dtType < 0
                            eIdx = eIdx + 1;
                        else
                            fIdx = fIdx + 1;
                        end
                        continue;
                    end
                    dt = (eegEventTimeSec(eIdx) - offset) - fnirsEventTimeSec(fIdx);
                    if abs(dt) <= alignToleranceSec
                        count = count + 1;
                        eIdx = eIdx + 1;
                        fIdx = fIdx + 1;
                    elseif dt < -alignToleranceSec
                        eIdx = eIdx + 1;
                    else
                        fIdx = fIdx + 1;
                    end
                end
                if count > bestMatch
                    bestMatch = count;
                    bestOffset = offset;
                end
            end
            offsetSec = bestOffset;

            eIdx = 1;
            fIdx = 1;
            while eIdx <= nEegEvents && fIdx <= nFnirsEvents
                if ~isfinite(eegEventTimeSec(eIdx))
                    eIdx = eIdx + 1;
                    continue;
                end
                if eegEventType(eIdx) ~= fnirsEventType(fIdx)
                    dtType = (eegEventTimeSec(eIdx) - offsetSec) - fnirsEventTimeSec(fIdx);
                    if dtType < 0
                        eIdx = eIdx + 1;
                    else
                        fIdx = fIdx + 1;
                    end
                    continue;
                end
                dt = (eegEventTimeSec(eIdx) - offsetSec) - fnirsEventTimeSec(fIdx);
                if abs(dt) <= alignToleranceSec
                    mapE(end+1,1) = eIdx; %#ok<AGROW>
                    mapF(end+1,1) = fIdx; %#ok<AGROW>
                    eIdx = eIdx + 1;
                    fIdx = fIdx + 1;
                elseif dt < -alignToleranceSec
                    eIdx = eIdx + 1;
                else
                    fIdx = fIdx + 1;
                end
            end
        else
            fprintf('[%s] WARNING: Not enough EEG event times; falling back to order.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            nMatch = min(nEegEvents, nFnirsEvents);
            mapE = (1:nMatch)';
            mapF = (1:nMatch)';
            if any(isfinite(eegEventTimeSec(1:nMatch)))
                offsetSec = median(eegEventTimeSec(1:nMatch) - fnirsEventTimeSec(1:nMatch), 'omitnan');
            else
                offsetSec = 0;
            end
        end
    end

    if ~isfinite(offsetSec)
        offsetSec = 0;
    end

    fnirsEventValidEEG = nan(nFnirsEvents, 1);
    for mi = 1:numel(mapE)
        eIdx = mapE(mi);
        fIdx = mapF(mi);
        fnirsEventValidEEG(fIdx) = double(eegEventValid(eIdx));
    end
    fnirsEventValidBoth = fnirsEventValid & (fnirsEventValidEEG == 1);

    matchedEvents = numel(mapE);
    matchedValid = sum(fnirsEventValidBoth);

    fprintf('[%s] EEG valid common=%d/%d oddball=%d/%d; fNIRS valid common=%d/%d oddball=%d/%d; matched=%d valid=%d\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        eegCommonValid, eegCommonCount, eegOddValid, eegOddCount, ...
        fnirsCommonValid, nCommonFnirs, fnirsOddValid, nOddFnirs, matchedEvents, matchedValid);

    alignRows(end+1, :) = {subName, eegCommonCount, eegCommonValid, eegOddCount, eegOddValid, ...
        nCommonFnirs, fnirsCommonValid, nOddFnirs, fnirsOddValid, matchedEvents, matchedValid, offsetSec};

    %% Load EEG continuous dataset and compute band power
    EEG_cont = pop_loadset('filename', eegContFile, 'filepath', eegSubDir);
    if EEG_cont.trials > 1
        fprintf('[%s] WARNING: EEG continuous data appears epoched for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end

    chanIdx = [];
    if isfield(EEG_cont, 'chanlocs') && ~isempty(EEG_cont.chanlocs)
        labels = {EEG_cont.chanlocs.labels};
        for li = 1:numel(eegChanLabels)
            matchIdx = find(strcmpi(labels, eegChanLabels{li}), 1);
            if ~isempty(matchIdx)
                chanIdx(end+1) = matchIdx; %#ok<AGROW>
            end
        end
    end
    if isempty(chanIdx)
        if eegChanFallbackIdx <= EEG_cont.nbchan
            chanIdx = eegChanFallbackIdx;
            fprintf('[%s] WARNING: EEG channel label not found; using index %d.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), chanIdx);
        else
            fprintf('[%s] WARNING: EEG channel not found and fallback invalid. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            continue;
        end
    end

    EEG_sel = pop_select(EEG_cont, 'channel', chanIdx);
    if isempty(EEG_sel.times)
        eegTimeSec = (0:EEG_sel.pnts-1)' / EEG_sel.srate;
    else
        eegTimeSec = double(EEG_sel.times(:)) / 1000;
    end

    numBands = size(bandDefs, 1);
    eegBandPower = cell(numBands, 1);
    for bi = 1:numBands
        band = bandDefs(bi, :);

        if exist('pop_eegfiltnew', 'file') == 2
            EEG_filt = pop_eegfiltnew(EEG_sel, 'locutoff', band(1), 'hicutoff', band(2), 'plotfreqz', 0);
        else
            if exist('eegfilt', 'file') ~= 2
                error('Neither pop_eegfiltnew nor eegfilt is available on the MATLAB path.');
            end
            EEG_filt = EEG_sel;
            data_filt = zeros(size(EEG_sel.data), 'like', EEG_sel.data);
            for ch = 1:size(EEG_sel.data, 1)
                data_filt(ch, :) = eegfilt(EEG_sel.data(ch, :), EEG_sel.srate, band(1), band(2));
            end
            EEG_filt.data = data_filt;
        end

        data = double(EEG_filt.data);
        analytic = hilbert(data') ;
        pow = abs(analytic).^2;

        if smoothWindowSec > 0
            winSamples = max(1, round(smoothWindowSec * EEG_sel.srate));
            pow = movmean(pow, winSamples, 1);
        end

        switch lower(eegCombineMode)
            case 'median'
                eegBandPower{bi} = median(pow, 2);
            otherwise
                eegBandPower{bi} = mean(pow, 2);
        end
    end

    %% fNIRS channel selection (HbO/HbR)
    ml = dc.measurementList;
    nMeasAll = numel(ml);
    if nMeasAll < 1
        fprintf('[%s] WARNING: fNIRS measurement list empty for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

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

    mlAct = true(nMeasAll,1);
    if isfield(S, 'mlActAuto') && ~isempty(S.mlActAuto)
        tmp = S.mlActAuto;
        if iscell(tmp)
            tmp = tmp{1};
        end
        tmp = double(tmp(:) > 0);
        if numel(tmp) == nMeasAll
            mlAct = tmp;
        end
    end

    roiNames = {};
    roiMasks = {};
    if useRoiFilter
        for ri = 1:numel(roiDefs)
            pairs = roiDefs(ri).pairs;
            mask = false(nMeasAll,1);
            for pi = 1:size(pairs,1)
                mask = mask | (src == pairs(pi,1) & det == pairs(pi,2));
            end
            roiNames{end+1} = roiDefs(ri).name; %#ok<AGROW>
            roiMasks{end+1} = mask; %#ok<AGROW>
        end
    else
        roiNames = {'All'};
        roiMasks = {true(nMeasAll,1)};
    end

    %% Resample EEG power to fNIRS time
    eegTimeAligned = eegTimeSec - offsetSec;
    fsFnirs = 1 / median(diff(tFnirs));
    if ~isfinite(fsFnirs) || fsFnirs <= 0
        fprintf('[%s] WARNING: Invalid fNIRS sampling rate for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    maxLagSamples = max(1, round(maxLagSec * fsFnirs));
    lagSecVec = (-maxLagSamples:maxLagSamples)' / fsFnirs;
    baseWinSamples = max(4, round(cohWindowSec * fsFnirs));
    baseNfft = max(cohNfft, baseWinSamples);

    eegBandPowerResamp = cell(numBands, 1);
    for bi = 1:numBands
        eegBandPowerResamp{bi} = interp1(eegTimeAligned, eegBandPower{bi}, tFnirs, 'linear', NaN);
    end

    %% Analysis masks
    baseMask = true(size(tFnirs));
    if ~isempty(tIncVec)
        baseMask = baseMask & (tIncVec > 0);
    end

    maskAll = false(size(tFnirs));
    maskCommon = false(size(tFnirs));
    maskOddball = false(size(tFnirs));

    if strcmpi(maskMode, 'matched_trials')
        for ei = 1:numel(fnirsEventTimeSec)
            if ~fnirsEventValidBoth(ei)
                continue;
            end
            winStart = fnirsEventTimeSec(ei) + analysisWindowSec(1);
            winStop = fnirsEventTimeSec(ei) + analysisWindowSec(2);
            idx = tFnirs >= winStart & tFnirs <= winStop;
            maskAll = maskAll | idx;
            if fnirsEventType(ei) == 1
                maskCommon = maskCommon | idx;
            elseif fnirsEventType(ei) == 2
                maskOddball = maskOddball | idx;
            end
        end
        maskAll = maskAll & baseMask;
        maskCommon = maskCommon & baseMask;
        maskOddball = maskOddball & baseMask;
    else
        maskAll = baseMask;
        maskCommon = baseMask;
        maskOddball = baseMask;
        fprintf('[%s] WARNING: maskMode=%s uses full time series; common/oddball analyses are identical.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), maskMode);
    end

    %% Cross-correlation and coherence
    chromList = {'HbO', 'HbR'};
    for ri = 1:numel(roiNames)
        roiMask = roiMasks{ri} & mlAct;
        if ~any(roiMask)
            fprintf('[%s] WARNING: ROI %s has no active channels for %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiNames{ri}, subName);
            continue;
        end
        for ci = 1:numel(chromList)
            chromMask = roiMask & (chrom == chromList{ci});
            if ~any(chromMask)
                fprintf('[%s] WARNING: ROI %s has no %s channels for %s.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiNames{ri}, chromList{ci}, subName);
                continue;
            end
            hbSeries = mean(yFnirs(:, chromMask), 2);

            for bi = 1:numBands
                eegSeries = eegBandPowerResamp{bi};
                condNames = {'All','Common','Oddball'};
                condMasks = {maskAll, maskCommon, maskOddball};
                condValid = false(1, numel(condNames));
                condSamples = zeros(1, numel(condNames));
                condPeakR = nan(1, numel(condNames));
                condPeakLag = nan(1, numel(condNames));
                condPeakIdx = nan(1, numel(condNames));
                condCohMean = nan(1, numel(condNames));
                condCohPeak = nan(1, numel(condNames));
                condCohPeakHz = nan(1, numel(condNames));
                condXc = cell(1, numel(condNames));
                condCoh = cell(1, numel(condNames));
                condF = cell(1, numel(condNames));

                for cj = 1:numel(condNames)
                    validMask = condMasks{cj} & isfinite(eegSeries) & isfinite(hbSeries);
                    condSamples(cj) = sum(validMask);
                    if condSamples(cj) < 10
                        fprintf('[%s] WARNING: Not enough samples for %s %s %s %s in %s.\n', ...
                            datestr(now, 'yyyy-mm-dd HH:MM:SS'), condNames{cj}, roiNames{ri}, chromList{ci}, bandNames{bi}, subName);
                        continue;
                    end

                    x = eegSeries(validMask);
                    y = hbSeries(validMask);

                    switch lower(normalizeMode)
                        case 'zscore'
                            x = (x - mean(x)) / (std(x) + eps);
                            y = (y - mean(y)) / (std(y) + eps);
                        case 'demean'
                            x = x - mean(x);
                            y = y - mean(y);
                        otherwise
                            % no normalization
                    end

                    xc = xcorr(x, y, maxLagSamples, 'coeff');
                    condXc{cj} = xc;
                    if strcmpi(xcorrPeakMode, 'abs')
                        [peakVal, idx] = max(abs(xc));
                        condPeakR(cj) = sign(xc(idx)) * peakVal;
                    else
                        [condPeakR(cj), idx] = max(xc);
                    end
                    condPeakLag(cj) = lagSecVec(idx);
                    condPeakIdx(cj) = idx;

                    winSamples = min(baseWinSamples, numel(x));
                    noverlap = min(winSamples-1, round(cohOverlapFrac * winSamples));
                    nfft = baseNfft;
                    [cxy, f] = mscohere(x, y, winSamples, noverlap, nfft, fsFnirs);
                    condCoh{cj} = cxy;
                    condF{cj} = f;

                    freqMask = f >= cohFreqRange(1) & f <= cohFreqRange(2);
                    if any(freqMask)
                        condCohMean(cj) = mean(cxy(freqMask));
                        [condCohPeak(cj), pIdx] = max(cxy(freqMask));
                        fSel = f(freqMask);
                        condCohPeakHz(cj) = fSel(pIdx);
                    else
                        condCohMean(cj) = NaN;
                        condCohPeak(cj) = NaN;
                        condCohPeakHz(cj) = NaN;
                    end

                    condValid(cj) = true;
                end

                outCondIdx = [1 2];
                outAnalysisIdx = [1 2];
                for oi = 1:numel(outCondIdx)
                    cj = outCondIdx(oi);
                    ai = outAnalysisIdx(oi);
                    if ~condValid(cj)
                        continue;
                    end
                    analysisName = analysisNames{ai};
                    analysisTag = analysisTags{ai};
                    analysisTitle = analysisTitles{ai};

                    summaryRows(end+1, :) = {subName, analysisName, roiNames{ri}, chromList{ci}, bandNames{bi}, ...
                        sprintf('%g-%g', bandDefs(bi,1), bandDefs(bi,2)), condSamples(cj), condPeakR(cj), condPeakLag(cj), ...
                        condCohMean(cj), condCohPeak(cj), condCohPeakHz(cj), offsetSec, matchedEvents, matchedValid};

                    if saveFullSpectra
                        outMat = fullfile(subOutDir, sprintf('%s_task-%s_step3_%s_%s_%s_%s.mat', ...
                            subName, taskName, analysisTag, roiNames{ri}, chromList{ci}, bandNames{bi}));
                        result = struct();
                        result.sub = subName;
                        result.analysis = analysisName;
                        result.roi = roiNames{ri};
                        result.chrom = chromList{ci};
                        result.band = bandNames{bi};
                        result.bandHz = bandDefs(bi,:);
                        result.offsetSec = offsetSec;
                        result.lagSec = lagSecVec;
                        result.xcorr = condXc{cj};
                        result.cohFreq = condF{cj};
                        result.coh = condCoh{cj};
                        result.samplesUsed = condSamples(cj);
                        save(outMat, 'result', '-v7.3');
                    end

                    if saveFigures && plotPerPairFigures && figCount < maxFiguresPerSubject
                        figCount = figCount + 1;
                        fig = figure('Visible', figureVisible);
                        subplot(2,1,1);
                        plot(lagSecVec, condXc{cj}, 'LineWidth', 1);
                        hold on;
                        plot(condPeakLag(cj), condXc{cj}(condPeakIdx(cj)), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
                        yl = ylim;
                        line([0 0], yl, 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
                        line([condPeakLag(cj) condPeakLag(cj)], yl, 'Color', [0.2 0.2 0.2], 'LineStyle', ':');
                        xlabel('Lag (s)');
                        ylabel('r');
                        title(sprintf('%s %s %s %s %s xcorr', subName, analysisTitle, roiNames{ri}, chromList{ci}, bandNames{bi}));
                        grid on;

                        subplot(2,1,2);
                        plot(condF{cj}, condCoh{cj}, 'LineWidth', 1);
                        hold on;
                        plot(condCohPeakHz(cj), condCohPeak(cj), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
                        xlabel('Hz');
                        ylabel('Coherence');
                        xlim([max(0, cohFreqRange(1)) cohFreqRange(2)]);
                        title(sprintf('%s coherence peak %.3f at %.3f Hz', analysisTitle, condCohPeak(cj), condCohPeakHz(cj)));
                        grid on;

                        baseName = sprintf('%s_%s_%s_%s_%s_xcorr_coh', subName, analysisTag, roiNames{ri}, chromList{ci}, bandNames{bi});
                        baseName = regexprep(baseName, '\s+', '');
                        for fi = 1:numel(figureFormats)
                            fmt = lower(figureFormats{fi});
                            outFile = fullfile(figSubDir, sprintf('%s.%s', baseName, fmt));
                            switch fmt
                                case 'fig'
                                    savefig(fig, outFile);
                                case 'png'
                                    print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
                                case 'jpg'
                                    print(fig, outFile, '-djpeg', sprintf('-r%d', figureDpi));
                                otherwise
                                    print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
                            end
                        end

                        if closeFigures
                            close(fig);
                        end
                    end
                end

                if condValid(2) && condValid(3)
                    if numel(condXc{2}) ~= numel(condXc{3}) || numel(condCoh{2}) ~= numel(condCoh{3})
                        fprintf('[%s] WARNING: Oddball-common length mismatch for %s %s %s in %s.\n', ...
                            datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiNames{ri}, chromList{ci}, bandNames{bi}, subName);
                    else
                        xcDiff = condXc{3} - condXc{2};
                        cohDiff = condCoh{3} - condCoh{2};
                        if strcmpi(xcorrPeakMode, 'abs')
                            [peakVal, idx] = max(abs(xcDiff));
                            peakR = sign(xcDiff(idx)) * peakVal;
                        else
                            [peakR, idx] = max(xcDiff);
                        end
                        peakLag = lagSecVec(idx);

                        f = condF{2};
                        freqMask = f >= cohFreqRange(1) & f <= cohFreqRange(2);
                        if any(freqMask)
                            cohMean = mean(cohDiff(freqMask));
                            [cohPeak, pIdx] = max(cohDiff(freqMask));
                            fSel = f(freqMask);
                            cohPeakHz = fSel(pIdx);
                        else
                            cohMean = NaN;
                            cohPeak = NaN;
                            cohPeakHz = NaN;
                        end

                        samplesUsed = min(condSamples(2), condSamples(3));
                        summaryRows(end+1, :) = {subName, analysisNames{3}, roiNames{ri}, chromList{ci}, bandNames{bi}, ...
                            sprintf('%g-%g', bandDefs(bi,1), bandDefs(bi,2)), samplesUsed, peakR, peakLag, ...
                            cohMean, cohPeak, cohPeakHz, offsetSec, matchedEvents, matchedValid};

                        if saveFullSpectra
                            outMat = fullfile(subOutDir, sprintf('%s_task-%s_step3_%s_%s_%s_%s.mat', ...
                                subName, taskName, analysisTags{3}, roiNames{ri}, chromList{ci}, bandNames{bi}));
                            result = struct();
                            result.sub = subName;
                            result.analysis = analysisNames{3};
                            result.roi = roiNames{ri};
                            result.chrom = chromList{ci};
                            result.band = bandNames{bi};
                            result.bandHz = bandDefs(bi,:);
                            result.offsetSec = offsetSec;
                            result.lagSec = lagSecVec;
                            result.xcorr = xcDiff;
                            result.cohFreq = f;
                            result.coh = cohDiff;
                            result.samplesUsed = samplesUsed;
                            save(outMat, 'result', '-v7.3');
                        end

                        if saveFigures && plotPerPairFigures && figCount < maxFiguresPerSubject
                            figCount = figCount + 1;
                            fig = figure('Visible', figureVisible);
                            subplot(2,1,1);
                            plot(lagSecVec, xcDiff, 'LineWidth', 1);
                            hold on;
                            plot(peakLag, xcDiff(idx), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
                            yl = ylim;
                            line([0 0], yl, 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
                            line([peakLag peakLag], yl, 'Color', [0.2 0.2 0.2], 'LineStyle', ':');
                            xlabel('Lag (s)');
                            ylabel('r');
                            title(sprintf('%s %s %s %s xcorr', subName, analysisTitles{3}, roiNames{ri}, chromList{ci}, bandNames{bi}));
                            grid on;

                            subplot(2,1,2);
                            plot(f, cohDiff, 'LineWidth', 1);
                            hold on;
                            plot(cohPeakHz, cohPeak, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
                            xlabel('Hz');
                            ylabel('Coherence');
                            xlim([max(0, cohFreqRange(1)) cohFreqRange(2)]);
                            title(sprintf('%s coherence peak %.3f at %.3f Hz', analysisTitles{3}, cohPeak, cohPeakHz));
                            grid on;

                            baseName = sprintf('%s_%s_%s_%s_%s_xcorr_coh', subName, analysisTags{3}, roiNames{ri}, chromList{ci}, bandNames{bi});
                            baseName = regexprep(baseName, '\s+', '');
                            for fi = 1:numel(figureFormats)
                                fmt = lower(figureFormats{fi});
                                outFile = fullfile(figSubDir, sprintf('%s.%s', baseName, fmt));
                                switch fmt
                                    case 'fig'
                                        savefig(fig, outFile);
                                    case 'png'
                                        print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
                                    case 'jpg'
                                        print(fig, outFile, '-djpeg', sprintf('-r%d', figureDpi));
                                    otherwise
                                        print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
                                end
                            end

                            if closeFigures
                                close(fig);
                            end
                        end
                    end
                else
                    fprintf('[%s] WARNING: Missing oddball/common samples for oddball-common in %s %s %s (%s).\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiNames{ri}, chromList{ci}, bandNames{bi}, subName);
                end
            end
        end
    end
end

%% Write summary tables
summaryTbl = [];
if ~isempty(summaryRows)
    summaryTbl = cell2table(summaryRows, 'VariableNames', summaryHeader);
end

if saveFigures && plotSummaryFigures && ~isempty(summaryTbl)
    roiNamePlot = summaryRoiName;
    roiList = unique(summaryTbl.ROI);
    if ~any(strcmpi(roiList, summaryRoiName)) && ~isempty(roiList)
        roiNamePlot = roiList{1};
        fprintf('[%s] WARNING: summaryRoiName=%s not found; using %s for plots.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), summaryRoiName, roiNamePlot);
    end

    for ai = 1:numel(analysisNames)
        analysisName = analysisNames{ai};
        analysisTag = analysisTags{ai};
        analysisTitle = analysisTitles{ai};
        roiMask = strcmpi(summaryTbl.ROI, roiNamePlot);
        analysisMask = strcmpi(summaryTbl.Analysis, analysisName);
        summaryRoiTbl = summaryTbl(roiMask & analysisMask, :);
        if isempty(summaryRoiTbl)
            fprintf('[%s] WARNING: No summary rows for %s (ROI=%s).\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName, roiNamePlot);
            continue;
        end

        peakLag = summaryRoiTbl.PeakLagSec;
        peakR = summaryRoiTbl.PeakR;
        if iscell(peakLag)
            peakLag = cell2mat(peakLag);
        end
        if iscell(peakR)
            peakR = cell2mat(peakR);
        end

        fig = figure('Visible', figureVisible);
        nBands = numel(bandNames);
        nChrom = numel(summaryChromList);
        for bi = 1:nBands
            for ci = 1:nChrom
                subplot(nBands, nChrom, (bi-1)*nChrom + ci);
                bandMask = strcmpi(summaryRoiTbl.Band, bandNames{bi});
                chromMask = strcmpi(summaryRoiTbl.Chrom, summaryChromList{ci});
                mask = bandMask & chromMask;
                if any(mask)
                    scatter(peakLag(mask), peakR(mask), 24, 'filled');
                else
                    scatter(0, 0, 4, 'filled');
                end
                xl = xlim;
                yl = ylim;
                line([0 0], yl, 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
                line(xl, [0 0], 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
                xlabel('Lag (s)');
                ylabel('r');
                title(sprintf('%s %s', bandNames{bi}, summaryChromList{ci}));
                grid on;
            end
        end
        if exist('sgtitle', 'file') == 2
            sgtitle(sprintf('Peak r vs lag (%s, ROI=%s)', analysisTitle, roiNamePlot));
        end

        baseName = sprintf('summary_peakr_lag_%s_roi-%s', analysisTag, roiNamePlot);
        baseName = regexprep(baseName, '\s+', '');
        for fi = 1:numel(figureFormats)
            fmt = lower(figureFormats{fi});
            outFile = fullfile(figRoot, sprintf('%s.%s', baseName, fmt));
            switch fmt
                case 'fig'
                    savefig(fig, outFile);
                case 'png'
                    print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
                case 'jpg'
                    print(fig, outFile, '-djpeg', sprintf('-r%d', figureDpi));
                otherwise
                    print(fig, outFile, '-dpng', sprintf('-r%d', figureDpi));
            end
        end

        if closeFigures
            close(fig);
        end
    end
end

summaryFile = fullfile(outRoot, 'STEP3_summary.csv');
fid = fopen(summaryFile, 'w');
if fid == -1
    error('Could not open summary file: %s', summaryFile);
end
fprintf(fid, '%s\n', strjoin(summaryHeader, ','));
for i = 1:size(summaryRows, 1)
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d\n', ...
        summaryRows{i,1}, summaryRows{i,2}, summaryRows{i,3}, summaryRows{i,4}, summaryRows{i,5}, ...
        summaryRows{i,6}, summaryRows{i,7}, summaryRows{i,8}, summaryRows{i,9}, summaryRows{i,10}, ...
        summaryRows{i,11}, summaryRows{i,12}, summaryRows{i,13}, summaryRows{i,14}, summaryRows{i,15});
end
fclose(fid);

alignFile = fullfile(outRoot, 'STEP3_alignment_summary.csv');
fid = fopen(alignFile, 'w');
if fid == -1
    error('Could not open alignment file: %s', alignFile);
end
fprintf(fid, '%s\n', strjoin(alignHeader, ','));
for i = 1:size(alignRows, 1)
    fprintf(fid, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.6f\n', ...
        alignRows{i,1}, alignRows{i,2}, alignRows{i,3}, alignRows{i,4}, alignRows{i,5}, ...
        alignRows{i,6}, alignRows{i,7}, alignRows{i,8}, alignRows{i,9}, ...
        alignRows{i,10}, alignRows{i,11}, alignRows{i,12});
end
fclose(fid);

fprintf('[%s] Done. Wrote %s and %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), summaryFile, alignFile);
