% STEP2_eeg_informed_fnirs.m
% EEG-informed fNIRS GLM for motor task using mu/beta ERD time courses.
% Assumptions & Notes:
% - EEG epochs are time-locked to block start markers (1=left, 2=right).
% - EEG and fNIRS block order is aligned within condition by onset order.
% - EEG ERD is computed from contralateral ROI (left hand -> right ROI, right hand -> left ROI).
% - Missing EEG blocks are excluded from fNIRS timepoints if excludeMissingEegBlocks=true.
% - fNIRS preproc MAT contains dc, stim, mlActAuto, tIncAuto, tIncAutoCh.

%% Config
close all; clearvars; clc;

taskName = 'motor';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset pattern priority (first match wins)
eegSetPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set', ...
    };

% fNIRS preproc MAT naming
fnirsPreprocPattern = '%s_task-%s_fnirs_preproc.mat';

% EEG condition mapping
cond_labels = {'left','right'};
cond_markers = [1 2];

% EEG ROI labels
roi_left_labels = {'C3','C1','C5','CP3','FC3'};
roi_right_labels = {'C4','C2','C6','FC4','CP4'};

% EEG ERD settings
baseline_ms = [-2000 0];
analysis_ms = [-2000 30000];
task_window_ms = [0 30000];
mu_band = [8 12];
beta_band = [13 30];
tfr_freqs = [4 40];
tfr_cycles = [3 0.5];
tfr_timesout = 200;
tfr_padratio = 2;

% fNIRS event settings
stimIndexLeft = 1;
stimIndexRight = 2;
stimDurationSec = 30;

% HRF settings (gamma)
hrfDurationSec = 30;
gammaShape = 6;
gammaScale = 1;

% GLM settings
glmMethod = 'ar-irls'; % 'ar-irls' | 'ols'
arIrlsOrderSec = 1;
driftModel = 'dct'; % 'none' | 'linear' | 'poly' | 'dct'
driftPolyOrder = 3;
dctCutoffSec = 256;

% EEG regressor options
includeTaskRegressors = true;
includeEegRegressors = true;
includeMu = true;
includeBeta = true;
orthogonalizeEegToTask = false;
zscoreEegRegressors = true;
excludeMissingEegBlocks = true;

% ROI filtering (same as fnirs/motor/motor_fnirs_processing/STEP4_GLM.m)
useRoiFilter = true;
saveRoiSummary = true;
useContraRoiOnly = true;
contraRoiForLeft = 'RightCentral';
contraRoiForRight = 'LeftCentral';
roiDefs = struct();
roiDefs(1).name = 'LeftCentral';
roiDefs(1).pairs = [ ...
    14 10; 10 10; 10 14; 14 14; 14 18; 18 18; 18 14];
roiDefs(2).name = 'RightCentral';
roiDefs(2).pairs = [ ...
    11 11; 15 11; 15 15; 11 15; 19 15; 19 19; 15 19];

% Output
outputDirName = 'STEP2_eeg_informed_fnirs';
saveFigures = true;
processOnlyFirst = false;

%% Setup
% scriptDir ='/Users/zhaifeifei/Desktop/eeg_fnirs/eeg_fnirs_processing/motor_eeg_fnirs_processing'

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', outputDirName);
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end

fprintf('[%s] Assumptions: markers {1=left,2=right}, block order matched within condition, contralateral ROI used.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
if excludeMissingEegBlocks
    fprintf('[%s] Missing EEG blocks are excluded from fNIRS timepoints.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
if exist('newtimef', 'file') ~= 2
    error('newtimef not found on MATLAB path (EEGLAB).');
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>

if strcmpi(glmMethod, 'ar-irls')
    if exist('ar_glm_final', 'file') ~= 2 || exist('robust_ar_fit', 'file') ~= 2
        error('AR-IRLS functions not found. Add Homer3 to path or set glmMethod="ols".');
    end
    if exist('robustfit', 'file') ~= 2 || exist('tcdf', 'file') ~= 2
        error('Statistics Toolbox required for AR-IRLS.');
    end
elseif strcmpi(glmMethod, 'ols')
    if exist('tcdf', 'file') ~= 2
        error('Statistics Toolbox required for OLS p-values (tcdf).');
    end
else
    error('Unknown glmMethod: %s', glmMethod);
end

%% Discover subjects
if isempty(subjList)
    eegSubs = dir(fullfile(eegRoot, 'sub-*'));
    fnirsSubs = dir(fullfile(fnirsRoot, 'sub-*'));
    eegNames = {eegSubs([eegSubs.isdir]).name};
    fnirsNames = {fnirsSubs([fnirsSubs.isdir]).name};
    subjList = intersect(eegNames, fnirsNames, 'stable');
end
if isempty(subjList)
    error('No subjects found under %s and %s', eegRoot, fnirsRoot);
end

summaryRows = cell(0, 10);

%% Main loop
for si = 1:numel(subjList)
    subName = subjList{si};
    fprintf('\n[%s] Subject %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        subName, si, numel(subjList));

    eegSubDir = fullfile(eegRoot, subName);
    fnirsSubDir = fullfile(fnirsRoot, subName);
    if exist(eegSubDir, 'dir') ~= 7 || exist(fnirsSubDir, 'dir') ~= 7
        fprintf('[%s] WARNING: Missing subject folder; skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end

    eegSetPath = '';
    eegSetFile = '';
    for pi = 1:numel(eegSetPatterns)
        candidate = sprintf(eegSetPatterns{pi}, subName, taskName);
        candidatePath = fullfile(eegSubDir, candidate);
        if exist(candidatePath, 'file') == 2
            eegSetPath = candidatePath;
            eegSetFile = candidate;
            break;
        end
    end
    if isempty(eegSetPath)
        fprintf('[%s] WARNING: EEG .set not found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    fnirsMat = fullfile(fnirsSubDir, sprintf(fnirsPreprocPattern, subName, taskName));
    if exist(fnirsMat, 'file') ~= 2
        tmp = dir(fullfile(fnirsSubDir, '*_fnirs_preproc.mat'));
        if isempty(tmp)
            fprintf('[%s] WARNING: fNIRS preproc MAT not found for %s. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
            continue;
        end
        fnirsMat = fullfile(tmp(1).folder, tmp(1).name);
    end

    subOutDir = fullfile(outRoot, subName);
    if exist(subOutDir, 'dir') ~= 7
        mkdir(subOutDir);
    end

    %% Load EEG
    EEG = pop_loadset('filename', eegSetFile, 'filepath', eegSubDir);
    if ~isfield(EEG, 'epoch') || EEG.trials < 1 || EEG.pnts < 2
        fprintf('[%s] WARNING: EEG data invalid for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % EEG ROI indices
    chan_labels = {EEG.chanlocs.labels};
    chan_labels = cellfun(@(x) upper(strtrim(x)), chan_labels, 'UniformOutput', false);
    roi_left_upper = cellfun(@upper, roi_left_labels, 'UniformOutput', false);
    roi_right_upper = cellfun(@upper, roi_right_labels, 'UniformOutput', false);
    left_idx = find(ismember(chan_labels, roi_left_upper));
    right_idx = find(ismember(chan_labels, roi_right_upper));
    if isempty(left_idx) || isempty(right_idx)
        fprintf('[%s] WARNING: EEG ROI channels missing for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % Epoch windows
    epoch_start_ms = EEG.xmin * 1000;
    epoch_end_ms = EEG.xmax * 1000;
    baseline_ms_use = [max(baseline_ms(1), epoch_start_ms), min(baseline_ms(2), min(0, epoch_end_ms))];
    analysis_ms_use = [max(analysis_ms(1), epoch_start_ms), min(analysis_ms(2), epoch_end_ms)];
    if baseline_ms_use(2) <= baseline_ms_use(1)
        baseline_ms_use = [epoch_start_ms, min(0, epoch_end_ms)];
        fprintf('[%s] WARNING: baseline window adjusted for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end
    if analysis_ms_use(2) <= analysis_ms_use(1)
        fprintf('[%s] WARNING: analysis window invalid for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    task_window_ms_use = [max(task_window_ms(1), analysis_ms_use(1)), min(task_window_ms(2), analysis_ms_use(2))];
    if task_window_ms_use(2) <= task_window_ms_use(1)
        task_window_ms_use = analysis_ms_use;
        fprintf('[%s] WARNING: task window adjusted to analysis window for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end

    % Good trial mask
    bad_trials = false(1, EEG.trials);
    if isfield(EEG, 'reject')
        rej_fields = fieldnames(EEG.reject);
        for f = 1:length(rej_fields)
            field_val = EEG.reject.(rej_fields{f});
            if isempty(field_val) || ~(islogical(field_val) || isnumeric(field_val))
                continue;
            end
            if isvector(field_val) && numel(field_val) == EEG.trials
                bad_trials = bad_trials | logical(field_val(:)');
            else
                field_size = size(field_val);
                if numel(field_size) == 2 && field_size(2) == EEG.trials
                    bad_trials = bad_trials | any(logical(field_val), 1);
                end
            end
        end
    end
    boundary_trials = false(1, EEG.trials);
    for ep = 1:EEG.trials
        if isfield(EEG.epoch(ep), 'eventtype')
            et = EEG.epoch(ep).eventtype;
            if iscell(et)
                boundary_trials(ep) = any(strcmpi(et, 'boundary'));
            elseif ischar(et)
                boundary_trials(ep) = strcmpi(et, 'boundary');
            end
        end
    end
    bad_trials = bad_trials | boundary_trials;
    good_trials = ~bad_trials;

    % Determine bin and event latency per epoch
    epoch_bin = nan(1, EEG.trials);
    epoch_event_latency = nan(1, EEG.trials);
    for ep = 1:EEG.trials
        if ~isfield(EEG.epoch(ep), 'event') || ~isfield(EEG.epoch(ep), 'eventlatency')
            continue;
        end
        event_idx = EEG.epoch(ep).event;
        event_lat = EEG.epoch(ep).eventlatency;

        if iscell(event_idx)
            event_idx_num = nan(1, numel(event_idx));
            for k = 1:numel(event_idx)
                if ischar(event_idx{k})
                    event_idx_num(k) = str2double(event_idx{k});
                else
                    event_idx_num(k) = double(event_idx{k});
                end
            end
            event_idx = event_idx_num;
        else
            if ischar(event_idx)
                event_idx = str2double(event_idx);
            else
                event_idx = double(event_idx);
            end
        end

        if iscell(event_lat)
            event_lat_num = nan(1, numel(event_lat));
            for k = 1:numel(event_lat)
                if ischar(event_lat{k})
                    event_lat_num(k) = str2double(event_lat{k});
                else
                    event_lat_num(k) = double(event_lat{k});
                end
            end
            event_lat = event_lat_num;
        else
            if ischar(event_lat)
                event_lat = str2double(event_lat);
            else
                event_lat = double(event_lat);
            end
        end

        if isempty(event_idx) || isempty(event_lat) || all(isnan(event_lat))
            continue;
        end

        [~, min_idx] = min(abs(event_lat));
        ev = event_idx(min_idx);
        if isnan(ev)
            continue;
        end

        bini = [];
        if isfield(EEG, 'EVENTLIST') && isfield(EEG.EVENTLIST, 'eventinfo')
            if ev >= 1 && ev <= length(EEG.EVENTLIST.eventinfo)
                bini = EEG.EVENTLIST.eventinfo(ev).bini;
            end
        end
        if ~isempty(bini)
            if iscell(bini)
                bini = bini{1};
            end
            if numel(bini) > 1
                bini = bini(1);
            end
            epoch_bin(ep) = double(bini);
        else
            if ev >= 1 && ev <= length(EEG.event)
                evtype = EEG.event(ev).type;
                if iscell(evtype)
                    evtype = evtype{1};
                end
                if ischar(evtype)
                    evtype_num = str2double(evtype);
                else
                    evtype_num = double(evtype);
                end
                if ~isnan(evtype_num)
                    epoch_bin(ep) = evtype_num;
                end
            end
        end

        if ev >= 1 && ev <= length(EEG.event)
            if isfield(EEG.event(ev), 'latency')
                epoch_event_latency(ep) = double(EEG.event(ev).latency);
            end
        end
    end

    % Block index per condition
    nCond = numel(cond_labels);
    cond_block_index = nan(nCond, EEG.trials);
    nBlocks_eeg = zeros(nCond, 1);
    for c = 1:nCond
        cond_trials = find(epoch_bin == cond_markers(c));
        if isempty(cond_trials)
            nBlocks_eeg(c) = 0;
            continue;
        end
        cond_lat = epoch_event_latency(cond_trials);
        valid_lat = isfinite(cond_lat);
        sorted_trials = [];
        if any(valid_lat)
            [~, sort_idx] = sort(cond_lat(valid_lat));
            sorted_trials = cond_trials(valid_lat);
            sorted_trials = sorted_trials(sort_idx);
        end
        if any(~valid_lat)
            sorted_trials = [sorted_trials, cond_trials(~valid_lat)];
        end
        cond_block_index(c, sorted_trials) = 1:numel(sorted_trials);
        nBlocks_eeg(c) = numel(sorted_trials);
    end

    % Compute ERD time courses per block
    eeg_times_ref = [];
    mu_sum = cell(nCond, 1);
    beta_sum = cell(nCond, 1);
    block_count = cell(nCond, 1);
    for c = 1:nCond
        mu_sum{c} = [];
        beta_sum{c} = [];
        block_count{c} = zeros(nBlocks_eeg(c), 1);
    end

    for tr = 1:EEG.trials
        if ~good_trials(tr)
            continue;
        end
        cond_idx = find(cond_markers == epoch_bin(tr), 1);
        if isempty(cond_idx)
            continue;
        end
        block_idx = cond_block_index(cond_idx, tr);
        if isnan(block_idx) || block_idx < 1
            continue;
        end
        if cond_idx == 1
            roi_idx = right_idx;
        else
            roi_idx = left_idx;
        end
        if isempty(roi_idx)
            continue;
        end

        roi_data = squeeze(mean(EEG.data(roi_idx, :, tr), 1));
        roi_data = double(roi_data(:));
        try
            [ersp, ~, ~, times, freqs] = newtimef(roi_data, EEG.pnts, ...
                [epoch_start_ms epoch_end_ms], EEG.srate, tfr_cycles, ...
                'baseline', baseline_ms_use, 'freqs', tfr_freqs, ...
                'timesout', tfr_timesout, 'padratio', tfr_padratio, ...
                'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
        catch ME
            fprintf('[%s] WARNING: newtimef failed for %s trial %d (%s). Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, tr, ME.message);
            continue;
        end

        time_mask = times >= analysis_ms_use(1) & times <= analysis_ms_use(2);
        if ~any(time_mask)
            continue;
        end
        mu_idx = freqs >= mu_band(1) & freqs <= mu_band(2);
        beta_idx = freqs >= beta_band(1) & freqs <= beta_band(2);
        if ~any(mu_idx) || ~any(beta_idx)
            continue;
        end
        mu_db = mean(ersp(mu_idx, time_mask), 1);
        beta_db = mean(ersp(beta_idx, time_mask), 1);
        mu_pct = (10 .^ (mu_db / 10) - 1) * 100;
        beta_pct = (10 .^ (beta_db / 10) - 1) * 100;
        times_sel = times(time_mask);

        if isempty(eeg_times_ref)
            eeg_times_ref = times_sel;
            nTimes = numel(eeg_times_ref);
            for c = 1:nCond
                mu_sum{c} = zeros(nBlocks_eeg(c), nTimes);
                beta_sum{c} = zeros(nBlocks_eeg(c), nTimes);
            end
        end

        if numel(times_sel) ~= numel(eeg_times_ref) || any(abs(times_sel - eeg_times_ref) > 1e-6)
            mu_pct = interp1(times_sel, mu_pct, eeg_times_ref, 'linear', 'extrap');
            beta_pct = interp1(times_sel, beta_pct, eeg_times_ref, 'linear', 'extrap');
        end

        if isempty(mu_sum{cond_idx})
            continue;
        end
        mu_sum{cond_idx}(block_idx, :) = mu_sum{cond_idx}(block_idx, :) + mu_pct;
        beta_sum{cond_idx}(block_idx, :) = beta_sum{cond_idx}(block_idx, :) + beta_pct;
        block_count{cond_idx}(block_idx) = block_count{cond_idx}(block_idx) + 1;
    end

    if isempty(eeg_times_ref)
        fprintf('[%s] WARNING: No EEG ERD time courses computed for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    mu_block = cell(nCond, 1);
    beta_block = cell(nCond, 1);
    blocks_good = zeros(nCond, 1);
    for c = 1:nCond
        if isempty(mu_sum{c})
            mu_block{c} = [];
            beta_block{c} = [];
            blocks_good(c) = 0;
            continue;
        end
        mu_block{c} = mu_sum{c};
        beta_block{c} = beta_sum{c};
        for b = 1:size(mu_block{c}, 1)
            if block_count{c}(b) > 0
                mu_block{c}(b, :) = mu_block{c}(b, :) ./ block_count{c}(b);
                beta_block{c}(b, :) = beta_block{c}(b, :) ./ block_count{c}(b);
            else
                mu_block{c}(b, :) = NaN;
                beta_block{c}(b, :) = NaN;
            end
        end
        blocks_good(c) = sum(block_count{c} > 0);
    end

    %% Load fNIRS
    S = load(fnirsMat, 'dc', 'stim', 'mlActAuto', 'tIncAuto', 'tIncAutoCh', 'qc');
    if ~isfield(S, 'dc') || ~isfield(S, 'stim')
        fprintf('[%s] WARNING: fNIRS preproc missing dc/stim for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    dc = S.dc;
    stim = S.stim;
    if iscell(stim)
        stim = [stim{:}];
    end

    t = double(dc.time);
    y = double(dc.dataTimeSeries);
    if size(y,1) ~= numel(t)
        fprintf('[%s] WARNING: fNIRS time length mismatch for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    nT = numel(t);

    if numel(stim) < max(stimIndexLeft, stimIndexRight)
        fprintf('[%s] WARNING: stim indices missing for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    stimLeft = stim(stimIndexLeft).data;
    stimRight = stim(stimIndexRight).data;
    if isempty(stimLeft) || isempty(stimRight)
        fprintf('[%s] WARNING: stim data empty for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    if size(stimLeft,2) == 1
        stimLeft = [stimLeft, repmat(stimDurationSec, size(stimLeft,1), 1), ones(size(stimLeft,1),1)];
    elseif size(stimLeft,2) == 2
        stimLeft = [stimLeft(:,1), stimLeft(:,2), ones(size(stimLeft,1),1)];
    else
        stimLeft = stimLeft(:,1:3);
    end
    if size(stimRight,2) == 1
        stimRight = [stimRight, repmat(stimDurationSec, size(stimRight,1), 1), ones(size(stimRight,1),1)];
    elseif size(stimRight,2) == 2
        stimRight = [stimRight(:,1), stimRight(:,2), ones(size(stimRight,1),1)];
    else
        stimRight = stimRight(:,1:3);
    end

    [~, sortIdxL] = sort(stimLeft(:,1));
    stimLeft = stimLeft(sortIdxL,:);
    [~, sortIdxR] = sort(stimRight(:,1));
    stimRight = stimRight(sortIdxR,:);

    nLeftFnirs = size(stimLeft, 1);
    nRightFnirs = size(stimRight, 1);

    %% HRF
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        fprintf('[%s] WARNING: Invalid fNIRS time vector for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end
    tHrf = (0:dt:hrfDurationSec)';
    hrf = (tHrf.^(gammaShape-1)) .* exp(-tHrf/gammaScale) ./ ...
        (gamma(gammaShape) * (gammaScale^gammaShape));
    if max(hrf) > 0
        hrf = hrf ./ max(hrf);
    end

    %% Task regressors
    uLeft = zeros(nT,1);
    for ei = 1:nLeftFnirs
        onset = stimLeft(ei,1);
        dur = stimLeft(ei,2);
        amp = stimLeft(ei,3);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        if ~isfinite(amp)
            amp = 1;
        end
        idx = t >= onset & t < (onset + dur);
        uLeft(idx) = uLeft(idx) + amp;
    end
    xLeft = conv(uLeft, hrf);
    xLeft = xLeft(1:nT);

    uRight = zeros(nT,1);
    for ei = 1:nRightFnirs
        onset = stimRight(ei,1);
        dur = stimRight(ei,2);
        amp = stimRight(ei,3);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        if ~isfinite(amp)
            amp = 1;
        end
        idx = t >= onset & t < (onset + dur);
        uRight(idx) = uRight(idx) + amp;
    end
    xRight = conv(uRight, hrf);
    xRight = xRight(1:nT);

    %% EEG-informed regressors (timecourse)
    eeg_times_sec = eeg_times_ref / 1000;
    analysis_window_sec = [analysis_ms_use(1) analysis_ms_use(2)] / 1000;
    task_mask = eeg_times_ref >= task_window_ms_use(1) & eeg_times_ref <= task_window_ms_use(2);
    if ~any(task_mask)
        task_mask = true(size(eeg_times_ref));
        fprintf('[%s] WARNING: task window outside ERD time axis for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
    end

    uEegMuLeft = zeros(nT,1);
    uEegBetaLeft = zeros(nT,1);
    uEegMuRight = zeros(nT,1);
    uEegBetaRight = zeros(nT,1);
    missingBlockMask = false(nT,1);

    alignRows = cell(0, 9);
    matchedLeft = 0;
    matchedRight = 0;

    for bi = 1:nLeftFnirs
        onset = stimLeft(bi,1);
        dur = stimLeft(bi,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        eegCount = 0;
        eegAvailable = false;
        if bi <= size(mu_block{1}, 1)
            eegCount = block_count{1}(bi);
            eegAvailable = eegCount > 0;
        end
        mu_mean = NaN;
        beta_mean = NaN;
        if eegAvailable
            mu_tc = mu_block{1}(bi, :);
            beta_tc = beta_block{1}(bi, :);
            idx = t >= onset + eeg_times_sec(1) & t <= onset + eeg_times_sec(end);
            if any(idx)
                rel = t(idx) - onset;
                mu_interp = interp1(eeg_times_sec, mu_tc, rel, 'linear', 'extrap');
                beta_interp = interp1(eeg_times_sec, beta_tc, rel, 'linear', 'extrap');
                uEegMuLeft(idx) = uEegMuLeft(idx) + mu_interp;
                uEegBetaLeft(idx) = uEegBetaLeft(idx) + beta_interp;
            end
            mu_mean = mean(mu_tc(task_mask), 'omitnan');
            beta_mean = mean(beta_tc(task_mask), 'omitnan');
            matchedLeft = matchedLeft + 1;
        else
            idxMissing = t >= onset + analysis_window_sec(1) & t <= onset + analysis_window_sec(2);
            missingBlockMask(idxMissing) = true;
        end
        alignRows(end+1,:) = {subName, cond_labels{1}, bi, onset, dur, eegAvailable, eegCount, mu_mean, beta_mean}; %#ok<AGROW>
    end

    for bi = 1:nRightFnirs
        onset = stimRight(bi,1);
        dur = stimRight(bi,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        eegCount = 0;
        eegAvailable = false;
        if bi <= size(mu_block{2}, 1)
            eegCount = block_count{2}(bi);
            eegAvailable = eegCount > 0;
        end
        mu_mean = NaN;
        beta_mean = NaN;
        if eegAvailable
            mu_tc = mu_block{2}(bi, :);
            beta_tc = beta_block{2}(bi, :);
            idx = t >= onset + eeg_times_sec(1) & t <= onset + eeg_times_sec(end);
            if any(idx)
                rel = t(idx) - onset;
                mu_interp = interp1(eeg_times_sec, mu_tc, rel, 'linear', 'extrap');
                beta_interp = interp1(eeg_times_sec, beta_tc, rel, 'linear', 'extrap');
                uEegMuRight(idx) = uEegMuRight(idx) + mu_interp;
                uEegBetaRight(idx) = uEegBetaRight(idx) + beta_interp;
            end
            mu_mean = mean(mu_tc(task_mask), 'omitnan');
            beta_mean = mean(beta_tc(task_mask), 'omitnan');
            matchedRight = matchedRight + 1;
        else
            idxMissing = t >= onset + analysis_window_sec(1) & t <= onset + analysis_window_sec(2);
            missingBlockMask(idxMissing) = true;
        end
        alignRows(end+1,:) = {subName, cond_labels{2}, bi, onset, dur, eegAvailable, eegCount, mu_mean, beta_mean}; %#ok<AGROW>
    end

    xEegMuLeft = conv(uEegMuLeft, hrf);
    xEegMuLeft = xEegMuLeft(1:nT);
    xEegMuRight = conv(uEegMuRight, hrf);
    xEegMuRight = xEegMuRight(1:nT);
    xEegBetaLeft = conv(uEegBetaLeft, hrf);
    xEegBetaLeft = xEegBetaLeft(1:nT);
    xEegBetaRight = conv(uEegBetaRight, hrf);
    xEegBetaRight = xEegBetaRight(1:nT);

    if zscoreEegRegressors
        reg = xEegMuLeft;
        mask = isfinite(reg) & reg ~= 0;
        if any(mask)
            mu_val = mean(reg(mask), 'omitnan');
            sd_val = std(reg(mask), 'omitnan');
            if ~isfinite(sd_val) || sd_val <= 0
                sd_val = 1;
            end
            xEegMuLeft = (reg - mu_val) / sd_val;
        end

        reg = xEegMuRight;
        mask = isfinite(reg) & reg ~= 0;
        if any(mask)
            mu_val = mean(reg(mask), 'omitnan');
            sd_val = std(reg(mask), 'omitnan');
            if ~isfinite(sd_val) || sd_val <= 0
                sd_val = 1;
            end
            xEegMuRight = (reg - mu_val) / sd_val;
        end

        reg = xEegBetaLeft;
        mask = isfinite(reg) & reg ~= 0;
        if any(mask)
            mu_val = mean(reg(mask), 'omitnan');
            sd_val = std(reg(mask), 'omitnan');
            if ~isfinite(sd_val) || sd_val <= 0
                sd_val = 1;
            end
            xEegBetaLeft = (reg - mu_val) / sd_val;
        end

        reg = xEegBetaRight;
        mask = isfinite(reg) & reg ~= 0;
        if any(mask)
            mu_val = mean(reg(mask), 'omitnan');
            sd_val = std(reg(mask), 'omitnan');
            if ~isfinite(sd_val) || sd_val <= 0
                sd_val = 1;
            end
            xEegBetaRight = (reg - mu_val) / sd_val;
        end
    end

    if orthogonalizeEegToTask && includeTaskRegressors
        Xtask = [xLeft xRight];
        if rank(Xtask) > 0
            xEegMuLeft = xEegMuLeft - Xtask * (Xtask \ xEegMuLeft);
            xEegMuRight = xEegMuRight - Xtask * (Xtask \ xEegMuRight);
            xEegBetaLeft = xEegBetaLeft - Xtask * (Xtask \ xEegBetaLeft);
            xEegBetaRight = xEegBetaRight - Xtask * (Xtask \ xEegBetaRight);
        end
    end

    %% Design matrix
    X = [];
    regNames = {};
    if includeTaskRegressors
        X = [X xLeft xRight];
        regNames = [regNames {'left','right'}];
    end
    if includeEegRegressors
        if includeMu
            X = [X xEegMuLeft xEegMuRight];
            regNames = [regNames {'eeg_mu_left','eeg_mu_right'}];
        end
        if includeBeta
            X = [X xEegBetaLeft xEegBetaRight];
            regNames = [regNames {'eeg_beta_left','eeg_beta_right'}];
        end
    end
    if isempty(X)
        fprintf('[%s] WARNING: No regressors selected for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    X = [X ones(nT,1)];
    regNames = [regNames {'intercept'}];

    switch lower(strtrim(driftModel))
        case 'none'
            % no drift terms
        case 'linear'
            tTrend = t - mean(t);
            X = [X tTrend];
            regNames = [regNames {'linear'}];
        case 'poly'
            if driftPolyOrder >= 1
                t0 = t - mean(t);
                tNorm = t0 / max(abs(t0));
                for po = 1:driftPolyOrder
                    X = [X tNorm.^po];
                    regNames = [regNames {sprintf('poly%d', po)}];
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
                        dctReg(:,k) = cos(pi * (n + 0.5) * k / nT);
                        regNames = [regNames {sprintf('dct%02d', k)}];
                    end
                    X = [X dctReg];
                end
            end
        otherwise
            error('Unknown driftModel: %s', driftModel);
    end

    %% ROI selection
    ml = dc.measurementList;
    nMeasAll = numel(ml);
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

    roiLabel = strings(nMeasAll,1);
    roiMask = true(nMeasAll,1);
    if useRoiFilter
        roiMask = false(nMeasAll,1);
        for ri = 1:numel(roiDefs)
            pairs = roiDefs(ri).pairs;
            pairMask = false(nMeasAll,1);
            for pi = 1:size(pairs,1)
                pairMask = pairMask | (src==pairs(pi,1) & det==pairs(pi,2));
            end
            roiLabel(pairMask) = roiDefs(ri).name;
            roiMask = roiMask | pairMask;
        end
        if ~any(roiMask)
            error('No channels matched ROI definitions. Check source/det pairs.');
        end
    else
        roiLabel(:) = "ALL";
    end

    if useContraRoiOnly && ~useRoiFilter
        warning('useContraRoiOnly requires useRoiFilter=true; disabling contra-only.');
        useContraRoiOnly = false;
    end
    contraLeftStr = string(contraRoiForLeft);
    contraRightStr = string(contraRoiForRight);

    %% Channel QC mask (mlActAuto)
    qcMaskAll = true(nMeasAll,1);
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
        if numel(mlActVec) == nMeasAll
            qcMaskAll = mlActVec;
        else
            warning('mlActAuto size mismatch (got %d, expected %d); QC mask disabled.', ...
                numel(mlActVec), nMeasAll);
        end
    else
        fprintf('[%s] WARNING: mlActAuto not found; using all channels.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    end

    glmMask = roiMask & qcMaskAll;
    if ~any(glmMask)
        error('No channels left after ROI+QC filtering.');
    end

    y = y(:, glmMask);
    src = src(glmMask);
    det = det(glmMask);
    chrom = chrom(glmMask);
    roiLabel = roiLabel(glmMask);
    nMeas = numel(src);

    %% Motion time masks (tIncAuto)
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
        if numel(tIncVec) ~= nT
            if numel(tIncVec) > nT
                tIncVec = tIncVec(1:nT);
                warning('tIncAuto longer than time vector; trimmed to %d.', nT);
            else
                padN = nT - numel(tIncVec);
                tIncVec = [tIncVec; ones(padN,1)];
                warning('tIncAuto shorter than time vector; padded to %d.', nT);
            end
        end
    end

    tIncCh = [];
    if isfield(S, 'tIncAutoCh') && ~isempty(S.tIncAutoCh)
        if iscell(S.tIncAutoCh)
            tIncCh = S.tIncAutoCh{1};
        else
            tIncCh = S.tIncAutoCh;
        end
    end
    if ~isempty(tIncCh)
        if size(tIncCh,1) ~= nT && size(tIncCh,2) == nT
            tIncCh = tIncCh';
        end
        if size(tIncCh,1) ~= nT
            warning('tIncAutoCh time dimension mismatch; ignoring channel-wise mask.');
            tIncCh = [];
        elseif size(tIncCh,2) ~= nMeasAll
            warning('tIncAutoCh channel dimension mismatch; ignoring channel-wise mask.');
            tIncCh = [];
        else
            tIncCh = double(tIncCh > 0);
            tIncCh = tIncCh(:, glmMask);
        end
    end

    %% GLM
    p = size(X,2);
    beta = nan(p, nMeas);
    seReg = nan(p, nMeas);
    tReg = nan(p, nMeas);
    pReg = nan(p, nMeas);
    dof = nan(1, nMeas);

    for mi = 1:nMeas
        tMask = true(nT,1);
        if ~isempty(tIncVec)
            tMask = tMask & (tIncVec > 0);
        end
        if ~isempty(tIncCh)
            tMask = tMask & (tIncCh(:,mi) > 0);
        end
        if excludeMissingEegBlocks && includeEegRegressors
            tMask = tMask & ~missingBlockMask;
        end
        tMask = tMask & all(isfinite(X), 2) & isfinite(y(:,mi));

        nKeep = sum(tMask);
        if nKeep <= p
            continue;
        end

        Xi = X(tMask, :);
        yi = y(tMask, mi);

        switch lower(glmMethod)
            case 'ar-irls'
                Pmax = max(1, round(arIrlsOrderSec / dt));
                [~, beta_i, tstat_i, pval_i, ~, CovB_i, dfe_i] = ar_glm_final(yi, Xi, Pmax);
                CovB_i = squeeze(CovB_i);
            case 'ols'
                dfe_i = size(Xi,1) - size(Xi,2);
                if dfe_i <= 0
                    continue;
                end
                beta_i = Xi \ yi;
                resid = yi - (Xi * beta_i);
                s2 = sum(resid.^2) / dfe_i;
                CovB_i = s2 * pinv(Xi' * Xi);
                se_i = sqrt(diag(CovB_i));
                tstat_i = beta_i ./ se_i;
                pval_i = 2 * tcdf(-abs(tstat_i), dfe_i);
            otherwise
                error('Unknown glmMethod: %s', glmMethod);
        end

        if isempty(CovB_i) || any(~isfinite(CovB_i(:)))
            continue;
        end

        beta(:, mi) = beta_i(:);
        seReg(:, mi) = sqrt(diag(CovB_i));
        tReg(:, mi) = tstat_i(:);
        pReg(:, mi) = pval_i(:);
        dof(mi) = dfe_i;
    end

    %% Output tables
    rowCell = {};
    for mi = 1:nMeas
        for ri = 1:numel(regNames)
            termName = regNames{ri};
            if useContraRoiOnly
                isLeftTerm = any(strcmp(termName, {'left','eeg_mu_left','eeg_beta_left'}));
                isRightTerm = any(strcmp(termName, {'right','eeg_mu_right','eeg_beta_right'}));
                if isLeftTerm && ~strcmp(roiLabel(mi), contraLeftStr)
                    continue;
                end
                if isRightTerm && ~strcmp(roiLabel(mi), contraRightStr)
                    continue;
                end
            end
            rowCell(end+1,:) = {subName, taskName, mi, src(mi), det(mi), char(chrom(mi)), char(roiLabel(mi)), ...
                termName, beta(ri,mi), seReg(ri,mi), tReg(ri,mi), pReg(ri,mi)}; %#ok<AGROW>
        end
    end

    T = cell2table(rowCell, 'VariableNames', ...
        {'subjId','task','measIndex','src','det','chromophore','roi','term','beta','SE','t','p'});

    outMat = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm.mat', subName, taskName));
    outCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm_table.csv', subName, taskName));
    save(outMat, 'beta', 'seReg', 'tReg', 'pReg', 'dof', 'regNames', 'X', ...
        'chrom', 'src', 'det', 'roiLabel', 'roiMask', 'glmMask', 'qcMaskAll', ...
        'eeg_times_ref', 'analysis_ms_use', 'task_window_ms_use', 'mu_band', 'beta_band', ...
        'block_count', 'missingBlockMask', '-v7.3');
    writetable(T, outCsv);

    %% ROI summary
    if saveRoiSummary
        if useRoiFilter
            roiList = {roiDefs.name};
        else
            roiList = {'ALL'};
        end
        chromOrder = ["HbO","HbR","HbT"];
        chromList = chromOrder(ismember(chromOrder, unique(chrom)));
        roiRows = {};
        for ri = 1:numel(roiList)
            roiName = roiList{ri};
            for ci = 1:numel(chromList)
                chromName = chromList(ci);
                chanIdx = (roiLabel == roiName) & (chrom == chromName);
                for ti2 = 1:numel(regNames)
                    termName = regNames{ti2};
                    if useContraRoiOnly
                        isLeftTerm = any(strcmp(termName, {'left','eeg_mu_left','eeg_beta_left'}));
                        isRightTerm = any(strcmp(termName, {'right','eeg_mu_right','eeg_beta_right'}));
                        if isLeftTerm && ~strcmp(roiName, contraRoiForLeft)
                            continue;
                        end
                        if isRightTerm && ~strcmp(roiName, contraRoiForRight)
                            continue;
                        end
                    end
                    regIdx = find(strcmp(regNames, termName), 1);
                    vals = beta(regIdx, chanIdx);
                    nChan = sum(chanIdx);
                    if nChan > 0
                        meanBeta = mean(vals, 'omitnan');
                        sdBeta = std(vals, 'omitnan');
                        seBeta = sdBeta / sqrt(nChan);
                    else
                        meanBeta = NaN;
                        sdBeta = NaN;
                        seBeta = NaN;
                    end
                    roiRows(end+1,:) = {subName, taskName, roiName, char(chromName), termName, nChan, meanBeta, sdBeta, seBeta}; %#ok<AGROW>
                end
            end
        end
        roiSummaryTable = cell2table(roiRows, 'VariableNames', ...
            {'subjId','task','roi','chromophore','term','nChannels','meanBeta','sdBeta','seBeta'});
        outRoiCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm_roi_table.csv', subName, taskName));
        writetable(roiSummaryTable, outRoiCsv);
        save(outMat, 'roiSummaryTable', '-append');
    end

    %% Alignment table
    alignTable = cell2table(alignRows, 'VariableNames', ...
        {'Subject','Condition','BlockIndex','FnirsOnsetSec','FnirsDurSec','EegBlockAvailable','EegGoodTrials','EegMuMean','EegBetaMean'});
    alignCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_fnirs_alignment.csv', subName, taskName));
    writetable(alignTable, alignCsv);

    %% Figures
    if saveFigures
        nPanels = 0;
        if includeTaskRegressors
            nPanels = nPanels + 1;
        end
        if includeMu
            nPanels = nPanels + 1;
        end
        if includeBeta
            nPanels = nPanels + 1;
        end
        if nPanels > 0
            figH = figure('Visible','off');
            panelIdx = 0;
            if includeTaskRegressors
                panelIdx = panelIdx + 1;
                subplot(nPanels,1,panelIdx);
                plot(t, xLeft, 'Color', [0.2 0.6 0.9]); hold on;
                plot(t, xRight, 'Color', [0.9 0.4 0.2]);
                title('Task regressors', 'Interpreter','none');
                legend({'Left','Right'}, 'Location','best');
            end
            if includeMu
                panelIdx = panelIdx + 1;
                subplot(nPanels,1,panelIdx);
                plot(t, xEegMuLeft, 'Color', [0.1 0.6 0.2]); hold on;
                plot(t, xEegMuRight, 'Color', [0.2 0.3 0.8]);
                title('EEG mu regressors (HRF-convolved)', 'Interpreter','none');
                legend({'Mu-left','Mu-right'}, 'Location','best');
            end
            if includeBeta
                panelIdx = panelIdx + 1;
                subplot(nPanels,1,panelIdx);
                plot(t, xEegBetaLeft, 'Color', [0.2 0.7 0.2]); hold on;
                plot(t, xEegBetaRight, 'Color', [0.6 0.2 0.8]);
                title('EEG beta regressors (HRF-convolved)', 'Interpreter','none');
                legend({'Beta-left','Beta-right'}, 'Location','best');
            end
            xlabel('Time (s)');
            figPath = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_regressors.png', subName, taskName));
            saveas(figH, figPath);
            close(figH);
        end
    end

    %% Summary row
    missingPct = 0;
    if excludeMissingEegBlocks && includeEegRegressors
        missingPct = 100 * mean(missingBlockMask);
    end
    summaryRows(end+1,:) = {subName, nBlocks_eeg(1), blocks_good(1), nBlocks_eeg(2), blocks_good(2), ...
        nLeftFnirs, nRightFnirs, matchedLeft, matchedRight, missingPct}; %#ok<AGROW>

    if processOnlyFirst
        break;
    end
end

%% Save summary
if ~isempty(summaryRows)
    summaryTable = cell2table(summaryRows, 'VariableNames', ...
        {'Subject','EEG_LeftBlocks','EEG_LeftGood','EEG_RightBlocks','EEG_RightGood', ...
        'FNIRS_LeftBlocks','FNIRS_RightBlocks','Matched_Left','Matched_Right','MissingTimePct'});
    summaryCsv = fullfile(outRoot, 'STEP2_eeg_informed_fnirs_summary.csv');
    writetable(summaryTable, summaryCsv);
end

fprintf('[%s] STEP2 EEG-informed fNIRS complete. Outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), outRoot);
