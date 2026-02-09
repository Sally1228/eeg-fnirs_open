% EEG-fNIRS dynamic coupling modeling (ARX / state-space / Kalman) for auditory oddball.
% Assumptions & Notes:
% - EEG bins: common=1 oddball=2 (update if your ERPLAB bin list differs).
% - fNIRS stim(1)=common stim(2)=oddball in the preproc MAT.
% - EEG and fNIRS alignment uses event timing when possible; fallback is order-based.
% - fNIRS trial validity uses tIncAuto (global); channel-wise tIncAutoCh is not used.
% - EEG input uses continuous band power (Hilbert) by default; P300 option uses oddball trials only.
% - Inputs/outputs are resampled to the fNIRS time grid; modeling is discrete-time at fNIRS sampling.
% - System Identification Toolbox is optional; if missing, ARX falls back to least squares.
% - Analyses include: All (common+oddball), Common-only, Oddball-Common (oddball minus common).
% - Oddball-Common uses event-locked oddball-minus-common windows defined by analysisWindowSec.

%% Config
close all; clearvars; clc;

taskName = 'audi';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset patterns (continuous and epoched)
eegContPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins.set', ...
    };
eegEpochPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set', ...
    };

% EEG bin mapping
eegCommonBin = 1;
eegOddballBin = 2;

% EEG channel selection
eegChanLabels = {'Pz'};
eegChanFallbackIdx = 53;
eegCombineMode = 'mean'; % 'mean' | 'median'

% EEG feature mode
eegFeatureMode = 'bandpower'; % 'bandpower' | 'p300'

% Band power settings
bandNames = {'Theta'};
bandDefs = [4 7];
powerMethod = 'hilbert'; % 'hilbert'
smoothEegSec = 0.0;

% P300 settings
p300WindowMs = [200 500];
baselineWindowMs = [-200 0];
p300Measure = 'peak'; % 'peak' | 'mean'
p300Polarity = 'positive'; % 'positive' | 'negative' | 'abs'
baselineMode = 'subtract'; % 'subtract' | 'none'
eventInputMode = 'box'; % 'box' | 'impulse'

% fNIRS stim settings
stimIndexCommon = 1;
stimIndexOddball = 2;
stimDurationSec = 0.5;
fnirsTrialKeepFrac = 0.8;
fnirsTrialPadSec = 0.0;

% Alignment settings
alignToleranceSec = 0.3;
minTimeMatchEvents = 5;

% Analysis conditions
analysisNames = {'All','Common','OddballCommon'};
analysisTags = {'all','common','oddballcommon'};
analysisWindowSec = [-2 20];
minSegmentSec = 6;
maskMode = 'matched_trials'; % 'matched_trials' | 'all' | 'tIncOnly'
applyTIncAuto = true;

% Resampling
resampleUniform = true;
uniformDtSec = [];

% Normalization/smoothing
normalizeMode = 'zscore'; % 'zscore' | 'demean' | 'none'
smoothFnirsSec = 0.0;

% Model settings
runArx = true;
runStateSpace = true;
runKalman = true;

arxNa = 2;
arxNb = 2;
delaySec = 2.5;
arxAllowLsFallback = true;

ssOrder = 2;

kalmanNa = 2;
kalmanNb = 2;
kalmanDelaySec = 2.5;
kalmanForgetFactor = 0.98;
kalmanQVar = 1e-4;
kalmanRVar = 1;
kalmanInitVar = 1e3;
kalmanBurnInSec = 5;
kalmanResetEachSegment = true;

irfDurationSec = 30;
irfPeakMode = 'abs'; % 'abs' | 'pos'

% fNIRS ROI settings
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

chromList = {'HbO'};

% Output
outputDirName = 'STEP4_dynamic_coupling';
savePerSubjectMat = true;
saveThetaTrace = false;
saveInputOutput = false;
processOnlyFirst = false;
saveExcel = true;
savePerSubjectExcel = true;

% Figure settings (summary)
saveFigures = true;
figureVisible = 'off'; % 'on' | 'off'
figureFormats = {'png'};
figureDpi = 150;
closeFigures = true;
maxSummaryFigures = 30;
plotMetric = 'Corr'; % 'Corr' | 'Gain' | 'RMSE' | 'IRFPeak' | 'IRFPeakLagSec' | 'DominantTauSec'
plotModels = {}; % empty = all
plotFeatures = {}; % empty = all
plotRoiList = {}; % empty = all
plotChromList = {}; % empty = all

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', outputDirName);
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end

if numel(analysisNames) ~= numel(analysisTags)
    error('analysisNames and analysisTags must have the same length.');
end
analysisLabel = strjoin(analysisNames, '/');

figRoot = fullfile(outRoot, 'figures');
if saveFigures && exist(figRoot, 'dir') ~= 7
    mkdir(figRoot);
end

fprintf('[%s] Assumptions: EEG bins common=%d oddball=%d; stim(1)=common stim(2)=oddball; tIncAuto trial QC; alignment by events; input=%s; analyses=%s; oddball-common uses event-locked oddball minus common windows.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), eegCommonBin, eegOddballBin, eegFeatureMode, analysisLabel);

if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %#ok<ASGLU>

if strcmpi(eegFeatureMode, 'bandpower') && strcmpi(powerMethod, 'hilbert') && exist('hilbert', 'file') ~= 2
    error('hilbert not found (Signal Processing Toolbox required).');
end
if strcmpi(eegFeatureMode, 'bandpower') && exist('pop_eegfiltnew', 'file') ~= 2 && exist('eegfilt', 'file') ~= 2
    error('Neither pop_eegfiltnew nor eegfilt is available on the MATLAB path.');
end

hasArx = exist('arx', 'file') == 2 && exist('iddata', 'file') == 2;
hasN4sid = exist('n4sid', 'file') == 2;
hasSsest = exist('ssest', 'file') == 2;

if runArx && ~hasArx && ~arxAllowLsFallback
    error('arx/iddata not found and LS fallback is disabled.');
end
if runStateSpace && ~(hasN4sid || hasSsest)
    fprintf('[%s] WARNING: n4sid/ssest not found; skipping state-space.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    runStateSpace = false;
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
summaryHeader = {'SubID','Analysis','ROI','Chrom','EEGFeature','BandHz','Model','Na','Nb','Nk', ...
    'SamplesUsed','SegmentsUsed','RMSE','Corr','Gain','IRFPeak','IRFPeakLagSec','DominantTauSec', ...
    'OffsetSec','MatchedEvents','MatchedValidEvents'};
summaryRows = cell(0, numel(summaryHeader));

alignHeader = {'SubID','EEGCommon','EEGCommonValid','EEGOdd','EEGOddValid', ...
    'FNIRSCommon','FNIRSCommonValid','FNIRSOdd','FNIRSOddValid', ...
    'MatchedEvents','MatchedValidEvents','OffsetSec'};
alignRows = cell(0, numel(alignHeader));

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

    %% Locate EEG datasets
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
    if strcmpi(eegFeatureMode, 'bandpower') && isempty(eegContPath)
        fprintf('[%s] WARNING: EEG continuous .set not found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

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
    eegEventTrialIdx = find(eegMask);

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

    dtRaw = median(diff(tFnirs));
    if resampleUniform
    dtTarget = uniformDtSec;
    if isempty(dtTarget)
        dtTarget = dtRaw;
    end
    tUniform = (tFnirs(1):dtTarget:tFnirs(end))';
    yFnirs = interp1(tFnirs, yFnirs, tUniform, 'linear', NaN);
    if ~isempty(tIncVec)
        tIncVec = interp1(tFnirs, tIncVec, tUniform, 'nearest', 0) > 0.5;
    end
    tFnirs = tUniform;
    end
    dt = median(diff(tFnirs));

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
            dtOffset = eegEventTimeSec - fnirsEventTimeSec;
            offsetSec = median(dtOffset(isfinite(dtOffset)));
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
                    dtMatch = (eegEventTimeSec(eIdx) - offset) - fnirsEventTimeSec(fIdx);
                    if abs(dtMatch) <= alignToleranceSec
                        count = count + 1;
                        eIdx = eIdx + 1;
                        fIdx = fIdx + 1;
                    elseif dtMatch < -alignToleranceSec
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
                dtMatch = (eegEventTimeSec(eIdx) - offsetSec) - fnirsEventTimeSec(fIdx);
                if abs(dtMatch) <= alignToleranceSec
                    mapE(end+1,1) = eIdx; %#ok<AGROW>
                    mapF(end+1,1) = fIdx; %#ok<AGROW>
                    eIdx = eIdx + 1;
                    fIdx = fIdx + 1;
                elseif dtMatch < -alignToleranceSec
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

    %% Build analysis masks
    maskAll = false(numel(tFnirs), 1);
    maskCommon = false(numel(tFnirs), 1);
    maskOddball = false(numel(tFnirs), 1);

    if strcmpi(maskMode, 'all') || strcmpi(maskMode, 'tinconly')
        maskAll = true(numel(tFnirs), 1);
        maskCommon = maskAll;
        maskOddball = maskAll;
    else
        for ei = 1:numel(fnirsEventTimeSec)
            if ~fnirsEventValidBoth(ei)
                continue;
            end
            winStart = fnirsEventTimeSec(ei) + analysisWindowSec(1);
            winStop = fnirsEventTimeSec(ei) + analysisWindowSec(2);
            idx = (tFnirs >= winStart & tFnirs <= winStop);
            maskAll = maskAll | idx;
            if fnirsEventType(ei) == 1
                maskCommon = maskCommon | idx;
            elseif fnirsEventType(ei) == 2
                maskOddball = maskOddball | idx;
            end
        end
    end
    if applyTIncAuto && ~isempty(tIncVec)
        maskAll = maskAll & tIncVec(:);
        maskCommon = maskCommon & tIncVec(:);
        maskOddball = maskOddball & tIncVec(:);
    end

    commonEventIdx = find(fnirsEventType == 1 & fnirsEventValidBoth);
    oddballEventIdx = find(fnirsEventType == 2 & fnirsEventValidBoth);

    windowSamples = max(1, round((analysisWindowSec(2) - analysisWindowSec(1)) / dt) + 1);
    windowDurationSec = (windowSamples - 1) * dt;
    tRel = (0:windowSamples-1)' * dt + analysisWindowSec(1);

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
            maskRoi = false(nMeasAll,1);
            for pi = 1:size(pairs,1)
                maskRoi = maskRoi | (src == pairs(pi,1) & det == pairs(pi,2));
            end
            roiNames{end+1} = roiDefs(ri).name; %#ok<AGROW>
            roiMasks{end+1} = maskRoi; %#ok<AGROW>
        end
    else
        roiNames = {'All'};
        roiMasks = {true(nMeasAll,1)};
    end

    %% EEG feature extraction
    uFullList = {};
    featureNames = {};
    featureBandHz = {};
    if strcmpi(eegFeatureMode, 'bandpower')
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
        for bi = 1:numBands
            band = bandDefs(bi, :);
            if exist('pop_eegfiltnew', 'file') == 2
                EEG_filt = pop_eegfiltnew(EEG_sel, 'locutoff', band(1), 'hicutoff', band(2), 'plotfreqz', 0);
            else
                EEG_filt = EEG_sel;
                data_filt = zeros(size(EEG_sel.data), 'like', EEG_sel.data);
                for ch = 1:size(EEG_sel.data, 1)
                    data_filt(ch, :) = eegfilt(EEG_sel.data(ch, :), EEG_sel.srate, band(1), band(2));
                end
                EEG_filt.data = data_filt;
            end

            data = double(EEG_filt.data);
            analytic = hilbert(data');
            pow = abs(analytic).^2;

            switch lower(eegCombineMode)
                case 'median'
                    eegPow = median(pow, 2);
                otherwise
                    eegPow = mean(pow, 2);
            end

            tEegAligned = eegTimeSec - offsetSec;
            uFull = interp1(tEegAligned, eegPow, tFnirs, 'linear', NaN);
            if smoothEegSec > 0
                winSamples = max(1, round(smoothEegSec / dt));
                uFull = movmean(uFull, winSamples);
            end

            uFullList{end+1} = uFull(:); %#ok<AGROW>
            featureNames{end+1} = bandNames{min(bi, numel(bandNames))}; %#ok<AGROW>
            featureBandHz{end+1} = sprintf('%g-%g', band(1), band(2)); %#ok<AGROW>
        end
    else
        chanIdx = [];
        if isfield(EEG_ep, 'chanlocs') && ~isempty(EEG_ep.chanlocs)
            labels = {EEG_ep.chanlocs.labels};
            for li = 1:numel(eegChanLabels)
                matchIdx = find(strcmpi(labels, eegChanLabels{li}), 1);
                if ~isempty(matchIdx)
                    chanIdx(end+1) = matchIdx; %#ok<AGROW>
                end
            end
        end
        if isempty(chanIdx)
            if eegChanFallbackIdx <= EEG_ep.nbchan
                chanIdx = eegChanFallbackIdx;
                fprintf('[%s] WARNING: EEG channel label not found; using index %d.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), chanIdx);
            else
                fprintf('[%s] WARNING: EEG channel not found and fallback invalid. Skipping.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
                continue;
            end
        end

        tMs = double(EEG_ep.times(:));
        winIdx = tMs >= p300WindowMs(1) & tMs <= p300WindowMs(2);
        baseIdx = tMs >= baselineWindowMs(1) & tMs <= baselineWindowMs(2);
        if ~any(winIdx)
            fprintf('[%s] WARNING: P300 window not found in EEG times for %s. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
            continue;
        end

        p300Amp = nan(EEG_ep.trials, 1);
        for tr = 1:EEG_ep.trials
            sig = double(EEG_ep.data(chanIdx, :, tr));
            if size(sig, 1) > 1
                sig = mean(sig, 1);
            end
            sig = sig(:)';
            if any(baseIdx)
                baselineVal = mean(sig(baseIdx));
            else
                baselineVal = 0;
            end
            if strcmpi(baselineMode, 'subtract')
                sig = sig - baselineVal;
            end
            seg = sig(winIdx);
            if isempty(seg)
                continue;
            end
            switch lower(p300Polarity)
                case 'negative'
                    val = min(seg);
                case 'abs'
                    val = max(abs(seg));
                otherwise
                    val = max(seg);
            end
            if strcmpi(p300Measure, 'mean')
                val = mean(seg);
            end
            p300Amp(tr) = val;
        end

        uFull = zeros(numel(tFnirs), 1);
        for mi = 1:numel(mapE)
            eIdx = mapE(mi);
            fIdx = mapF(mi);
            if eegEventType(eIdx) ~= 2
                continue;
            end
            if ~(eegEventValid(eIdx) && fnirsEventValid(fIdx))
                continue;
            end
            trialIdx = eegEventTrialIdx(eIdx);
            if trialIdx < 1 || trialIdx > numel(p300Amp)
                continue;
            end
            amp = p300Amp(trialIdx);
            if ~isfinite(amp)
                continue;
            end
            onset = fnirsEventTimeSec(fIdx);
            [~, idx] = min(abs(tFnirs - onset));
            if strcmpi(eventInputMode, 'box')
                durSamples = max(1, round(stimDurationSec / dt));
                idxEnd = min(numel(tFnirs), idx + durSamples - 1);
                uFull(idx:idxEnd) = uFull(idx:idxEnd) + amp;
            else
                uFull(idx) = uFull(idx) + amp;
            end
        end

        uFullList{1} = uFull(:);
        featureNames{1} = 'P300';
        featureBandHz{1} = 'P300';
    end

    if isempty(uFullList)
        fprintf('[%s] WARNING: No EEG features extracted for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    %% Loop over ROI, chrom, and feature
    modelResults = struct('analysis', {}, 'roi', {}, 'chrom', {}, 'feature', {}, 'bandHz', {}, ...
        'modelType', {}, 'A', {}, 'B', {}, 'C', {}, 'D', {}, 'nk', {}, ...
        'rmse', {}, 'corr', {}, 'gain', {}, 'irf', {}, 'irfTimeSec', {}, ...
        'irfPeak', {}, 'irfPeakLagSec', {}, 'dominantTauSec', {}, ...
        'samplesUsed', {}, 'segmentsUsed', {}, 'thetaMean', {}, 'thetaStd', {}, 'thetaTrace', {});

    for ri = 1:numel(roiNames)
        roiName = roiNames{ri};
        roiMask = roiMasks{ri};
        for ci = 1:numel(chromList)
            targetChrom = chromList{ci};
            chromMask = roiMask & mlAct & (chrom == string(targetChrom));
            if ~any(chromMask)
                fprintf('[%s] WARNING: No channels for ROI %s %s in %s. Skipping.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), roiName, targetChrom, subName);
                continue;
            end

            yFull = mean(yFnirs(:, chromMask), 2);
            if smoothFnirsSec > 0
                winSamples = max(1, round(smoothFnirsSec / dt));
                yFull = movmean(yFull, winSamples);
            end

            for fi = 1:numel(uFullList)
                uFull = uFullList{fi};

                if numel(uFull) ~= numel(yFull)
                    fprintf('[%s] WARNING: Length mismatch u/y for %s. Skipping.\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
                    continue;
                end

                for ai = 1:numel(analysisTags)
                    analysisTag = analysisTags{ai};
                    analysisName = analysisNames{ai};
                    useOddballCommon = strcmpi(analysisTag, 'oddballcommon');

                    if strcmpi(analysisTag, 'common')
                        analysisMask = maskCommon;
                    elseif strcmpi(analysisTag, 'oddball')
                        analysisMask = maskOddball;
                    else
                        analysisMask = maskAll;
                    end

                    goodMask = analysisMask & isfinite(uFull) & isfinite(yFull);
                    if sum(goodMask) < 10
                        fprintf('[%s] WARNING: Not enough valid samples for %s %s %s (%s). Skipping.\n', ...
                            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, targetChrom, analysisName);
                        continue;
                    end

                    uNorm = uFull;
                    yNorm = yFull;
                    if strcmpi(normalizeMode, 'zscore')
                        uMu = mean(uFull(goodMask));
                        uSd = std(uFull(goodMask));
                        yMu = mean(yFull(goodMask));
                        ySd = std(yFull(goodMask));
                        if uSd == 0 || ySd == 0
                            fprintf('[%s] WARNING: Zero std in normalization for %s (%s). Skipping.\n', ...
                                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, analysisName);
                            continue;
                        end
                        uNorm = (uFull - uMu) / uSd;
                        yNorm = (yFull - yMu) / ySd;
                    elseif strcmpi(normalizeMode, 'demean')
                        uMu = mean(uFull(goodMask));
                        yMu = mean(yFull(goodMask));
                        uNorm = uFull - uMu;
                        yNorm = yFull - yMu;
                    end

                    uSegs = {};
                    ySegs = {};
                    if useOddballCommon
                        if isempty(commonEventIdx) || isempty(oddballEventIdx)
                            fprintf('[%s] WARNING: Missing common/oddball events for %s %s %s (%s). Skipping.\n', ...
                                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, targetChrom, analysisName);
                            continue;
                        end
                        if windowDurationSec < minSegmentSec
                            fprintf('[%s] WARNING: Window shorter than minSegmentSec for %s %s %s (%s). Skipping.\n', ...
                                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, targetChrom, analysisName);
                            continue;
                        end

                        % Build oddball-minus-common segments using event-locked windows.
                        commonMatU = nan(windowSamples, numel(commonEventIdx));
                        commonMatY = nan(windowSamples, numel(commonEventIdx));
                        commonUsed = 0;
                        for ce = 1:numel(commonEventIdx)
                            ei = commonEventIdx(ce);
                            tAbs = fnirsEventTimeSec(ei) + tRel;
                            uSeg = interp1(tFnirs, uNorm, tAbs, 'linear', NaN);
                            ySeg = interp1(tFnirs, yNorm, tAbs, 'linear', NaN);
                            if applyTIncAuto && ~isempty(tIncVec)
                                incSeg = interp1(tFnirs, double(tIncVec), tAbs, 'nearest', 0);
                                uSeg(incSeg < 0.5) = NaN;
                                ySeg(incSeg < 0.5) = NaN;
                            end
                            if any(~isfinite(uSeg)) || any(~isfinite(ySeg))
                                continue;
                            end
                            commonUsed = commonUsed + 1;
                            commonMatU(:, commonUsed) = uSeg;
                            commonMatY(:, commonUsed) = ySeg;
                        end
                        if commonUsed < 1
                            fprintf('[%s] WARNING: No valid common segments for %s %s %s (%s). Skipping.\n', ...
                                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, targetChrom, analysisName);
                            continue;
                        end
                        commonMatU = commonMatU(:, 1:commonUsed);
                        commonMatY = commonMatY(:, 1:commonUsed);
                        uCommonMean = mean(commonMatU, 2);
                        yCommonMean = mean(commonMatY, 2);

                        for oe = 1:numel(oddballEventIdx)
                            ei = oddballEventIdx(oe);
                            tAbs = fnirsEventTimeSec(ei) + tRel;
                            uSeg = interp1(tFnirs, uNorm, tAbs, 'linear', NaN);
                            ySeg = interp1(tFnirs, yNorm, tAbs, 'linear', NaN);
                            if applyTIncAuto && ~isempty(tIncVec)
                                incSeg = interp1(tFnirs, double(tIncVec), tAbs, 'nearest', 0);
                                uSeg(incSeg < 0.5) = NaN;
                                ySeg(incSeg < 0.5) = NaN;
                            end
                            if any(~isfinite(uSeg)) || any(~isfinite(ySeg))
                                continue;
                            end
                            uDiff = uSeg - uCommonMean;
                            yDiff = ySeg - yCommonMean;
                            if any(~isfinite(uDiff)) || any(~isfinite(yDiff))
                                continue;
                            end
                            uSegs{end+1} = uDiff; %#ok<AGROW>
                            ySegs{end+1} = yDiff; %#ok<AGROW>
                        end
                    else
                        goodMask = analysisMask & isfinite(uNorm) & isfinite(yNorm);
                        edges = diff([0; goodMask; 0]);
                        segStarts = find(edges == 1);
                        segEnds = find(edges == -1) - 1;
                        for si2 = 1:numel(segStarts)
                            sIdx = segStarts(si2);
                            eIdx = segEnds(si2);
                            if (eIdx - sIdx + 1) * dt < minSegmentSec
                                continue;
                            end
                            uSegs{end+1} = uNorm(sIdx:eIdx); %#ok<AGROW>
                            ySegs{end+1} = yNorm(sIdx:eIdx); %#ok<AGROW>
                        end
                    end

                    if isempty(uSegs)
                        fprintf('[%s] WARNING: No valid segments for %s %s %s (%s). Skipping.\n', ...
                            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, roiName, targetChrom, analysisName);
                        continue;
                    end

                    samplesUsed = sum(cellfun(@numel, ySegs));
                    segmentsUsed = numel(ySegs);

                    %% ARX model
                    if runArx
                        nkArx = max(1, round(delaySec / dt));
                        maxLagArx = max(arxNa, arxNb + nkArx - 1);
                        if samplesUsed <= maxLagArx + 2
                            fprintf('[%s] WARNING: Too few samples for ARX in %s (%s).\n', ...
                                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, analysisName);
                        else
                            A = [];
                            B = [];
                            modelLabel = 'ARX';
                            if hasArx
                                data = iddata(ySegs, uSegs, dt);
                                arxModel = arx(data, [arxNa arxNb nkArx]);
                                A = arxModel.A(:).';
                                B = arxModel.B(:).';
                            elseif arxAllowLsFallback
                                modelLabel = 'ARX-LS';
                                nParams = arxNa + arxNb;
                                XtX = zeros(nParams);
                                Xty = zeros(nParams,1);
                                totalRows = 0;
                                for sj = 1:numel(ySegs)
                                    ySeg = ySegs{sj};
                                    uSeg = uSegs{sj};
                                    N = numel(ySeg);
                                    if N <= maxLagArx + 1
                                        continue;
                                    end
                                    nRows = N - maxLagArx;
                                    X = zeros(nRows, nParams);
                                    yVec = ySeg(maxLagArx+1:end);
                                    for t = 1:nRows
                                        idx = t + maxLagArx;
                                        phiY = -ySeg(idx-1:-1:idx-arxNa);
                                        phiU = uSeg(idx-nkArx:-1:idx-nkArx-arxNb+1);
                                        X(t,:) = [phiY(:)' phiU(:)'];
                                    end
                                    XtX = XtX + X' * X;
                                    Xty = Xty + X' * yVec;
                                    totalRows = totalRows + nRows;
                                end
                                if totalRows >= nParams
                                    theta = XtX \ Xty;
                                else
                                    theta = pinv(XtX) * Xty;
                                end
                                A = [1, theta(1:arxNa)'];
                                B = [zeros(1, nkArx-1), theta(arxNa+1:end)'];
                            end

                        if ~isempty(A) && ~isempty(B)
                            irfSamples = max(1, round(irfDurationSec / dt));
                            uImp = zeros(irfSamples, 1);
                            uImp(1) = 1;
                            irf = filter(B, A, uImp);
                            if strcmpi(irfPeakMode, 'pos')
                                [irfPeak, peakIdx] = max(irf);
                            else
                                [~, peakIdx] = max(abs(irf));
                                irfPeak = irf(peakIdx);
                            end
                            irfPeakLagSec = (peakIdx - 1) * dt;
                            irfTimeSec = (0:irfSamples-1)' * dt;

                            yAll = [];
                            yHatAll = [];
                            for sj = 1:numel(ySegs)
                                ySeg = ySegs{sj};
                                uSeg = uSegs{sj};
                                yHat = filter(B, A, uSeg);
                                startIdx = maxLagArx + 1;
                                if numel(ySeg) > startIdx
                                    yAll = [yAll; ySeg(startIdx:end)]; %#ok<AGROW>
                                    yHatAll = [yHatAll; yHat(startIdx:end)]; %#ok<AGROW>
                                end
                            end
                            if isempty(yAll)
                                rmse = NaN;
                                rVal = NaN;
                            else
                                rmse = sqrt(mean((yAll - yHatAll).^2));
                                R = corrcoef(yAll, yHatAll);
                                if size(R,1) > 1
                                    rVal = R(1,2);
                                else
                                    rVal = NaN;
                                end
                            end

                            sumA = sum(A);
                            if sumA ~= 0
                                gain = sum(B) / sumA;
                            else
                                gain = NaN;
                            end
                            poles = roots(A);
                            if all(abs(poles) < 1)
                                tau = -dt ./ log(abs(poles));
                                dominantTauSec = max(tau(isfinite(tau)));
                            else
                                dominantTauSec = NaN;
                            end

                            modelResults(end+1) = struct('analysis', analysisName, 'roi', roiName, 'chrom', targetChrom, ...
                                'feature', featureNames{fi}, 'bandHz', featureBandHz{fi}, ...
                                'modelType', modelLabel, 'A', A, 'B', B, 'C', [], 'D', [], 'nk', nkArx, ...
                                'rmse', rmse, 'corr', rVal, 'gain', gain, 'irf', irf, 'irfTimeSec', irfTimeSec, ...
                                'irfPeak', irfPeak, 'irfPeakLagSec', irfPeakLagSec, 'dominantTauSec', dominantTauSec, ...
                                'samplesUsed', samplesUsed, 'segmentsUsed', segmentsUsed, ...
                                'thetaMean', [], 'thetaStd', [], 'thetaTrace', []);

                            summaryRows(end+1, :) = {subName, analysisName, roiName, targetChrom, featureNames{fi}, featureBandHz{fi}, ...
                                modelLabel, arxNa, arxNb, nkArx, samplesUsed, segmentsUsed, rmse, rVal, gain, ...
                                irfPeak, irfPeakLagSec, dominantTauSec, offsetSec, matchedEvents, matchedValid};
                        end
                    end
                end

                %% State-space model
                if runStateSpace
                    if samplesUsed > ssOrder + 2
                        if hasN4sid
                            data = iddata(ySegs, uSegs, dt);
                            ssModel = n4sid(data, ssOrder, 'Ts', dt);
                        else
                            data = iddata(ySegs, uSegs, dt);
                            ssModel = ssest(data, ssOrder);
                        end

                        A = ssModel.A;
                        B = ssModel.B;
                        C = ssModel.C;
                        D = ssModel.D;

                        irfSamples = max(1, round(irfDurationSec / dt));
                        uImp = zeros(irfSamples, 1);
                        uImp(1) = 1;
                        x = zeros(size(A,1), 1);
                        irf = zeros(irfSamples, 1);
                        for t = 1:irfSamples
                            irf(t) = C * x + D * uImp(t);
                            x = A * x + B * uImp(t);
                        end
                        if strcmpi(irfPeakMode, 'pos')
                            [irfPeak, peakIdx] = max(irf);
                        else
                            [~, peakIdx] = max(abs(irf));
                            irfPeak = irf(peakIdx);
                        end
                        irfPeakLagSec = (peakIdx - 1) * dt;
                        irfTimeSec = (0:irfSamples-1)' * dt;

                        yAll = [];
                        yHatAll = [];
                        for sj = 1:numel(ySegs)
                            ySeg = ySegs{sj};
                            uSeg = uSegs{sj};
                            x = zeros(size(A,1), 1);
                            yHat = zeros(numel(ySeg), 1);
                            for t = 1:numel(ySeg)
                                yHat(t) = C * x + D * uSeg(t);
                                x = A * x + B * uSeg(t);
                            end
                            yAll = [yAll; ySeg(:)]; %#ok<AGROW>
                            yHatAll = [yHatAll; yHat(:)]; %#ok<AGROW>
                        end
                        rmse = sqrt(mean((yAll - yHatAll).^2));
                        R = corrcoef(yAll, yHatAll);
                        if size(R,1) > 1
                            rVal = R(1,2);
                        else
                            rVal = NaN;
                        end

                        if all(abs(eig(A)) < 1)
                            gain = C * ((eye(size(A)) - A) \ B) + D;
                            tau = -dt ./ log(abs(eig(A)));
                            dominantTauSec = max(tau(isfinite(tau)));
                        else
                            gain = NaN;
                            dominantTauSec = NaN;
                        end

                        modelResults(end+1) = struct('analysis', analysisName, 'roi', roiName, 'chrom', targetChrom, ...
                            'feature', featureNames{fi}, 'bandHz', featureBandHz{fi}, ...
                            'modelType', 'SS', 'A', A, 'B', B, 'C', C, 'D', D, 'nk', 0, ...
                            'rmse', rmse, 'corr', rVal, 'gain', gain, 'irf', irf, 'irfTimeSec', irfTimeSec, ...
                            'irfPeak', irfPeak, 'irfPeakLagSec', irfPeakLagSec, 'dominantTauSec', dominantTauSec, ...
                            'samplesUsed', samplesUsed, 'segmentsUsed', segmentsUsed, ...
                            'thetaMean', [], 'thetaStd', [], 'thetaTrace', []);

                        summaryRows(end+1, :) = {subName, analysisName, roiName, targetChrom, featureNames{fi}, featureBandHz{fi}, ...
                            'SS', ssOrder, NaN, 0, samplesUsed, segmentsUsed, rmse, rVal, gain, ...
                            irfPeak, irfPeakLagSec, dominantTauSec, offsetSec, matchedEvents, matchedValid};
                    end
                end

                %% Kalman time-varying ARX
                if runKalman
                    nkKal = max(1, round(kalmanDelaySec / dt));
                    maxLagKal = max(kalmanNa, kalmanNb + nkKal - 1);
                    if samplesUsed > maxLagKal + 2
                        nParams = kalmanNa + kalmanNb;
                        theta = zeros(nParams, 1);
                        P = eye(nParams) * kalmanInitVar;
                        burnInSamples = max(0, round(kalmanBurnInSec / dt));
                        thetaSum = zeros(nParams, 1);
                        thetaSumSq = zeros(nParams, 1);
                        thetaCount = 0;
                        yAll = [];
                        yHatAll = [];
                        thetaTrace = [];

                        for sj = 1:numel(ySegs)
                            ySeg = ySegs{sj};
                            uSeg = uSegs{sj};
                            if kalmanResetEachSegment
                                theta = zeros(nParams, 1);
                                P = eye(nParams) * kalmanInitVar;
                            end
                            for t = (maxLagKal+1):numel(ySeg)
                                phiY = -ySeg(t-1:-1:t-kalmanNa);
                                phiU = uSeg(t-nkKal:-1:t-nkKal-kalmanNb+1);
                                phi = [phiY(:); phiU(:)];
                                if any(~isfinite(phi))
                                    continue;
                                end
                                P = P / kalmanForgetFactor + kalmanQVar * eye(nParams);
                                Sval = phi' * P * phi + kalmanRVar;
                                K = (P * phi) / Sval;
                                yHat = phi' * theta;
                                err = ySeg(t) - yHat;
                                theta = theta + K * err;
                                P = P - K * phi' * P;
                                if t > maxLagKal + burnInSamples
                                    thetaSum = thetaSum + theta;
                                    thetaSumSq = thetaSumSq + theta.^2;
                                    thetaCount = thetaCount + 1;
                                end
                                yAll = [yAll; ySeg(t)]; %#ok<AGROW>
                                yHatAll = [yHatAll; yHat]; %#ok<AGROW>
                                if saveThetaTrace
                                    thetaTrace = [thetaTrace; theta']; %#ok<AGROW>
                                end
                            end
                        end

                        if thetaCount > 0
                            thetaMean = thetaSum / thetaCount;
                            thetaStd = sqrt(max(thetaSumSq / thetaCount - thetaMean.^2, 0));
                        else
                            thetaMean = theta;
                            thetaStd = NaN(size(theta));
                        end
                        A = [1, thetaMean(1:kalmanNa)'];
                        B = [zeros(1, nkKal-1), thetaMean(kalmanNa+1:end)'];

                        irfSamples = max(1, round(irfDurationSec / dt));
                        uImp = zeros(irfSamples, 1);
                        uImp(1) = 1;
                        irf = filter(B, A, uImp);
                        if strcmpi(irfPeakMode, 'pos')
                            [irfPeak, peakIdx] = max(irf);
                        else
                            [~, peakIdx] = max(abs(irf));
                            irfPeak = irf(peakIdx);
                        end
                        irfPeakLagSec = (peakIdx - 1) * dt;
                        irfTimeSec = (0:irfSamples-1)' * dt;

                        rmse = sqrt(mean((yAll - yHatAll).^2));
                        R = corrcoef(yAll, yHatAll);
                        if size(R,1) > 1
                            rVal = R(1,2);
                        else
                            rVal = NaN;
                        end

                        sumA = sum(A);
                        if sumA ~= 0
                            gain = sum(B) / sumA;
                        else
                            gain = NaN;
                        end
                        poles = roots(A);
                        if all(abs(poles) < 1)
                            tau = -dt ./ log(abs(poles));
                            dominantTauSec = max(tau(isfinite(tau)));
                        else
                            dominantTauSec = NaN;
                        end

                        thetaTraceOut = [];
                        if saveThetaTrace
                            thetaTraceOut = thetaTrace;
                        end

                        modelResults(end+1) = struct('analysis', analysisName, 'roi', roiName, 'chrom', targetChrom, ...
                            'feature', featureNames{fi}, 'bandHz', featureBandHz{fi}, ...
                            'modelType', 'Kalman', 'A', A, 'B', B, 'C', [], 'D', [], 'nk', nkKal, ...
                            'rmse', rmse, 'corr', rVal, 'gain', gain, 'irf', irf, 'irfTimeSec', irfTimeSec, ...
                            'irfPeak', irfPeak, 'irfPeakLagSec', irfPeakLagSec, 'dominantTauSec', dominantTauSec, ...
                            'samplesUsed', samplesUsed, 'segmentsUsed', segmentsUsed, ...
                            'thetaMean', thetaMean, 'thetaStd', thetaStd, 'thetaTrace', thetaTraceOut);

                        summaryRows(end+1, :) = {subName, analysisName, roiName, targetChrom, featureNames{fi}, featureBandHz{fi}, ...
                            'Kalman', kalmanNa, kalmanNb, nkKal, samplesUsed, segmentsUsed, rmse, rVal, gain, ...
                            irfPeak, irfPeakLagSec, dominantTauSec, offsetSec, matchedEvents, matchedValid};
                    end
                end
            end
        end
    end
    end

    %% Save per-subject results
    if savePerSubjectMat
        outMat = fullfile(subOutDir, sprintf('%s_task-%s_step4_dynamic_coupling.mat', subName, taskName));
        results = struct();
        results.subName = subName;
        results.taskName = taskName;
        results.offsetSec = offsetSec;
        results.matchedEvents = matchedEvents;
        results.matchedValid = matchedValid;
        results.modelResults = modelResults;
        if saveInputOutput
            results.uFullList = uFullList;
            results.featureNames = featureNames;
            results.featureBandHz = featureBandHz;
        end
        save(outMat, 'results');
    end

    if savePerSubjectExcel
        subRows = cell(0, numel(summaryHeader));
        for mr = 1:numel(modelResults)
            res = modelResults(mr);
            naVal = NaN;
            nbVal = NaN;
            if strcmpi(res.modelType, 'SS')
                naVal = ssOrder;
                nbVal = NaN;
            elseif strcmpi(res.modelType, 'Kalman')
                naVal = kalmanNa;
                nbVal = kalmanNb;
            else
                naVal = arxNa;
                nbVal = arxNb;
            end
            subRows(end+1, :) = {subName, res.analysis, res.roi, res.chrom, res.feature, res.bandHz, ...
                res.modelType, naVal, nbVal, res.nk, res.samplesUsed, res.segmentsUsed, res.rmse, res.corr, ...
                res.gain, res.irfPeak, res.irfPeakLagSec, res.dominantTauSec, offsetSec, matchedEvents, matchedValid}; %#ok<AGROW>
        end
        subTable = cell2table(subRows, 'VariableNames', summaryHeader);
        outXlsx = fullfile(subOutDir, sprintf('%s_task-%s_step4_dynamic_coupling.xlsx', subName, taskName));
        writetable(subTable, outXlsx, 'Sheet', 'Summary');
    end

    if processOnlyFirst
        break;
    end
end

%% Save summary tables
summaryTable = [];
if ~isempty(summaryRows)
    summaryTable = cell2table(summaryRows, 'VariableNames', summaryHeader);
    outCsv = fullfile(outRoot, 'STEP4_summary.csv');
    writetable(summaryTable, outCsv);
    if saveExcel
        outXlsx = fullfile(outRoot, 'STEP4_summary.xlsx');
        writetable(summaryTable, outXlsx);
    end
end

if ~isempty(alignRows)
    alignTable = cell2table(alignRows, 'VariableNames', alignHeader);
    outAlign = fullfile(outRoot, 'STEP4_alignment_summary.csv');
    writetable(alignTable, outAlign);
    if saveExcel
        outAlignXlsx = fullfile(outRoot, 'STEP4_alignment_summary.xlsx');
        writetable(alignTable, outAlignXlsx);
    end
end

%% Plot summary figures
if saveFigures && ~isempty(summaryTable)
    metricField = '';
    switch lower(plotMetric)
        case 'corr'
            metricField = 'Corr';
        case 'gain'
            metricField = 'Gain';
        case 'rmse'
            metricField = 'RMSE';
        case 'irfpeak'
            metricField = 'IRFPeak';
        case 'irfpeaklagsec'
            metricField = 'IRFPeakLagSec';
        case 'dominanttausec'
            metricField = 'DominantTauSec';
    end
    if isempty(metricField) || ~ismember(metricField, summaryTable.Properties.VariableNames)
        fprintf('[%s] WARNING: Plot metric %s not found in summary; skipping figures.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), plotMetric);
    else
        metricVals = summaryTable.(metricField);
        if iscell(metricVals)
            metricVals = cell2mat(metricVals);
        end

        modelList = unique(summaryTable.Model, 'stable');
        featureList = unique(summaryTable.EEGFeature, 'stable');
        roiList = unique(summaryTable.ROI, 'stable');
        chromListPlot = unique(summaryTable.Chrom, 'stable');

        if ~isempty(plotModels)
            keep = false(size(modelList));
            for i = 1:numel(modelList)
                keep(i) = any(strcmpi(modelList{i}, plotModels));
            end
            modelList = modelList(keep);
        end
        if ~isempty(plotFeatures)
            keep = false(size(featureList));
            for i = 1:numel(featureList)
                keep(i) = any(strcmpi(featureList{i}, plotFeatures));
            end
            featureList = featureList(keep);
        end
        if ~isempty(plotRoiList)
            keep = false(size(roiList));
            for i = 1:numel(roiList)
                keep(i) = any(strcmpi(roiList{i}, plotRoiList));
            end
            roiList = roiList(keep);
        end
        if ~isempty(plotChromList)
            keep = false(size(chromListPlot));
            for i = 1:numel(chromListPlot)
                keep(i) = any(strcmpi(chromListPlot{i}, plotChromList));
            end
            chromListPlot = chromListPlot(keep);
        end

        figCount = 0;
        for mi = 1:numel(modelList)
            for fi = 1:numel(featureList)
                figCount = figCount + 1;
                if figCount > maxSummaryFigures
                    break;
                end
                fig = figure('Visible', figureVisible);
                nRoi = max(1, numel(roiList));
                nChrom = max(1, numel(chromListPlot));
                for ri = 1:nRoi
                    for ci = 1:nChrom
                        subplot(nRoi, nChrom, (ri-1)*nChrom + ci);
                        baseMask = strcmpi(summaryTable.Model, modelList{mi}) & ...
                            strcmpi(summaryTable.EEGFeature, featureList{fi}) & ...
                            strcmpi(summaryTable.ROI, roiList{ri}) & ...
                            strcmpi(summaryTable.Chrom, chromListPlot{ci});
                        means = nan(1, numel(analysisNames));
                        sems = nan(1, numel(analysisNames));
                        valsByCond = cell(1, numel(analysisNames));
                        for ai = 1:numel(analysisNames)
                            condMask = strcmpi(summaryTable.Analysis, analysisNames{ai});
                            vals = metricVals(baseMask & condMask);
                            vals = vals(isfinite(vals));
                            valsByCond{ai} = vals;
                            if ~isempty(vals)
                                means(ai) = mean(vals);
                                sems(ai) = std(vals) / sqrt(numel(vals));
                            end
                        end
                        hold on;
                        bar(1:numel(analysisNames), means, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.4 0.4 0.4]);
                        errorbar(1:numel(analysisNames), means, sems, 'k.', 'LineWidth', 1);
                        for ai = 1:numel(analysisNames)
                            vals = valsByCond{ai};
                            if isempty(vals)
                                continue;
                            end
                            if numel(vals) > 1
                                jitter = linspace(-0.07, 0.07, numel(vals));
                            else
                                jitter = 0;
                            end
                            scatter(ai + jitter(:), vals(:), 18, 'filled');
                        end
                        set(gca, 'XTick', 1:numel(analysisNames), 'XTickLabel', analysisNames);
                        if exist('xtickangle', 'file') == 2
                            xtickangle(30);
                        end
                        ylabel(metricField);
                        title(sprintf('%s %s', roiList{ri}, chromListPlot{ci}));
                        grid on;
                    end
                end
                if exist('sgtitle', 'file') == 2
                    sgtitle(sprintf('%s %s %s', metricField, modelList{mi}, featureList{fi}));
                end

                baseName = sprintf('summary_%s_%s_%s', lower(metricField), modelList{mi}, featureList{fi});
                baseName = regexprep(baseName, '\s+', '');
                for fi2 = 1:numel(figureFormats)
                    fmt = lower(figureFormats{fi2});
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
            if figCount > maxSummaryFigures
                break;
            end
        end
    end
end
