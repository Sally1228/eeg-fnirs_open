% STEP2_eeg_informed_fnirs.m
% EEG-informed fNIRS GLM for cate task using Delta/Theta/Alpha/Beta ERD time courses.
% Assumptions & Notes:
% - EEG epochs are time-locked to marker=1 task blocks.
% - EEG and fNIRS block order is aligned by onset order.
% - EEG ERD is computed from ROI-specific channels (LeftFrontal/LeftTemporal/RightFrontal/RightTemporal).
% - Delta/Theta/Alpha/Beta are modeled in separate GLMs (one band per model).
% - Missing EEG blocks can be excluded from fNIRS timepoints if excludeMissingEegBlocks=true.
% - fNIRS preproc MAT contains dc, stim, mlActAuto, tIncAuto, tIncAutoCh.

%% Config
close all; clearvars; clc;

taskName = 'cate';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset pattern priority (first match wins)
eegSetPatterns = {'%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'};

% fNIRS preproc MAT naming
fnirsPreprocPattern = '%s_task-%s_fnirs_preproc.mat';

% Single condition mapping
cond_marker = 1;

% EEG ROI labels
roi_names = {'LeftFrontal','LeftTemporal','RightFrontal','RightTemporal'};
roi_labels = {
    {'F1','F3','F5','AF3','FC1','FC3','FC5','C1','C3'}, ...
    {'FT7','T7','TP7','CP5','C5'}, ...
    {'F2','F4','F6','AF4','FC2','FC4','FC6','C2','C4'}, ...
    {'FT8','T8','TP8','CP6','C6'} ...
};

% EEG ERD settings
baseline_ms = [-2000 0];
analysis_ms = [-2000 60000];
task_window_ms = [0 60000];
band_defs = {
    'Delta', [1 4];
    'Theta', [4 8];
    'Alpha', [8 13];
    'Beta',  [13 30]
};
band_names = band_defs(:,1);
nBand = size(band_defs, 1);
tfr_freqs = [1 40];
tfr_cycles = [3 0.5];
tfr_timesout = 200;
tfr_padratio = 2;

% fNIRS event settings
stimIndexTask = 1;
stimDurationSec = 60;

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
orthogonalizeEegToTask = false;
zscoreEegRegressors = true;
excludeMissingEegBlocks = true;

% ROI filtering (same as fnirs/cate/cate_fnirs_processing/STEP4_GLM.m)
useRoiFilter = true;
saveRoiSummary = true;
useRoiSpecificOnly = true;
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

% Output
outputDirName = 'STEP2_eeg_informed_fnirs';
saveFigures = true;
processOnlyFirst = false;

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', outputDirName);
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end

fprintf('[%s] Assumptions: marker=1 task blocks, EEG/fNIRS block order matched by onset.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Assumptions: EEG ROI=Left/Right Frontal/Temporal; EEG band ERD used as regressors.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Bands: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), strjoin(band_names', ', '));
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

summaryRows = cell(0, 8);

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
    roi_idx_all = cell(1, numel(roi_names));
    for r = 1:numel(roi_names)
        roi_upper = cellfun(@upper, roi_labels{r}, 'UniformOutput', false);
        roi_idx_all{r} = find(ismember(chan_labels, roi_upper));
        if isempty(roi_idx_all{r})
            fprintf('[%s] WARNING: EEG ROI %s channels missing for %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), roi_names{r}, subName);
        end
    end
    if all(cellfun(@isempty, roi_idx_all))
        fprintf('[%s] WARNING: No EEG ROI channels found for %s. Skipping.\n', ...
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

    % Block index for task marker
    cond_trials = find(epoch_bin == cond_marker);
    if isempty(cond_trials)
        fprintf('[%s] WARNING: no marker=%d epochs for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), cond_marker, subName);
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

    cond_block_index = nan(1, EEG.trials);
    cond_block_index(sorted_trials) = 1:numel(sorted_trials);
    nBlocks_eeg = numel(sorted_trials);

    % Compute ERD time courses per block and ROI (all bands)
    eeg_times_ref = [];
    nRoi = numel(roi_names);
    band_sum = cell(nRoi, nBand);
    band_count = cell(nRoi, nBand);
    for r = 1:nRoi
        for bi = 1:nBand
            band_sum{r,bi} = [];
            band_count{r,bi} = zeros(nBlocks_eeg, 1);
        end
    end

    for tr = 1:EEG.trials
        if ~good_trials(tr)
            continue;
        end
        if epoch_bin(tr) ~= cond_marker
            continue;
        end
        block_idx = cond_block_index(tr);
        if isnan(block_idx) || block_idx < 1
            continue;
        end

        for r = 1:nRoi
            roi_idx = roi_idx_all{r};
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
                fprintf('[%s] WARNING: newtimef failed for %s trial %d ROI %s (%s). Skipping.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, tr, roi_names{r}, ME.message);
                continue;
            end

            time_mask = times >= analysis_ms_use(1) & times <= analysis_ms_use(2);
            if ~any(time_mask)
                continue;
            end
            times_sel = times(time_mask);
            band_tc = nan(nBand, numel(times_sel));
            for bi = 1:nBand
                band_range = band_defs{bi,2};
                freq_mask = freqs >= band_range(1) & freqs <= band_range(2);
                if any(freq_mask)
                    band_db = mean(ersp(freq_mask, time_mask), 1);
                    band_tc(bi, :) = (10 .^ (band_db / 10) - 1) * 100;
                end
            end

            if isempty(eeg_times_ref)
                eeg_times_ref = times_sel;
                nTimes = numel(eeg_times_ref);
                for rr = 1:nRoi
                    for bi = 1:nBand
                        band_sum{rr,bi} = zeros(nBlocks_eeg, nTimes);
                    end
                end
            end

            if numel(times_sel) ~= numel(eeg_times_ref) || any(abs(times_sel - eeg_times_ref) > 1e-6)
                for bi = 1:nBand
                    band_tc(bi, :) = interp1(times_sel, band_tc(bi, :), eeg_times_ref, 'linear', 'extrap');
                end
            end

            for bi = 1:nBand
                if isempty(band_sum{r,bi})
                    continue;
                end
                if all(isfinite(band_tc(bi, :)))
                    band_sum{r,bi}(block_idx, :) = band_sum{r,bi}(block_idx, :) + band_tc(bi, :);
                    band_count{r,bi}(block_idx) = band_count{r,bi}(block_idx) + 1;
                end
            end
        end
    end

    if isempty(eeg_times_ref)
        fprintf('[%s] WARNING: No EEG ERD time courses computed for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    band_block = cell(nRoi, nBand);
    blocks_good = zeros(nRoi, nBand);
    for r = 1:nRoi
        for bi = 1:nBand
            if isempty(band_sum{r,bi})
                band_block{r,bi} = [];
                blocks_good(r,bi) = 0;
                continue;
            end
            band_block{r,bi} = band_sum{r,bi};
            for b = 1:size(band_block{r,bi}, 1)
                if band_count{r,bi}(b) > 0
                    band_block{r,bi}(b, :) = band_block{r,bi}(b, :) ./ band_count{r,bi}(b);
                else
                    band_block{r,bi}(b, :) = NaN;
                end
            end
            blocks_good(r,bi) = sum(band_count{r,bi} > 0);
        end
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

    if numel(stim) < stimIndexTask
        fprintf('[%s] WARNING: stim index missing for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    stimTask = stim(stimIndexTask).data;
    if isempty(stimTask)
        fprintf('[%s] WARNING: stim data empty for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    if size(stimTask,2) == 1
        stimTask = [stimTask, repmat(stimDurationSec, size(stimTask,1), 1), ones(size(stimTask,1),1)];
    elseif size(stimTask,2) == 2
        stimTask = [stimTask(:,1), stimTask(:,2), ones(size(stimTask,1),1)];
    else
        stimTask = stimTask(:,1:3);
    end

    [~, sortIdx] = sort(stimTask(:,1));
    stimTask = stimTask(sortIdx,:);
    nTaskFnirs = size(stimTask, 1);

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

    %% Task regressor
    uTask = zeros(nT,1);
    for ei = 1:nTaskFnirs
        onset = stimTask(ei,1);
        dur = stimTask(ei,2);
        amp = stimTask(ei,3);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        if ~isfinite(amp)
            amp = 1;
        end
        idx = t >= onset & t < (onset + dur);
        uTask(idx) = uTask(idx) + amp;
    end
    xTask = conv(uTask, hrf);
    xTask = xTask(1:nT);

    %% EEG-informed regressors (timecourse)
    eeg_times_sec = eeg_times_ref / 1000;
    analysis_window_sec = [analysis_ms_use(1) analysis_ms_use(2)] / 1000;
    task_mask = eeg_times_ref >= task_window_ms_use(1) & eeg_times_ref <= task_window_ms_use(2);
    if ~any(task_mask)
        task_mask = true(size(eeg_times_ref));
        fprintf('[%s] WARNING: task window outside ERD time axis for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
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

    %% Band-specific GLM loop
    for bi = 1:nBand
        band_name = band_defs{bi,1};
        band_tag = lower(band_name);
        band_range = band_defs{bi,2};

        uEegBand = zeros(nT, nRoi);
        missingBlockMask = false(nT,1);

        alignRows = cell(0, 9);
        matchedBlocks = zeros(nRoi, 1);

        for blk = 1:nTaskFnirs
            onset = stimTask(blk,1);
            dur = stimTask(blk,2);
            if ~isfinite(dur) || dur <= 0
                dur = stimDurationSec;
            end
            for r = 1:nRoi
                eegCount = 0;
                eegAvailable = false;
                if blk <= size(band_block{r,bi}, 1)
                    eegCount = band_count{r,bi}(blk);
                    eegAvailable = eegCount > 0;
                end
                band_mean = NaN;
                if eegAvailable
                    band_tc = band_block{r,bi}(blk, :);
                    idx = t >= onset + eeg_times_sec(1) & t <= onset + eeg_times_sec(end);
                    if any(idx)
                        rel = t(idx) - onset;
                        band_interp = interp1(eeg_times_sec, band_tc, rel, 'linear', 'extrap');
                        uEegBand(idx, r) = uEegBand(idx, r) + band_interp;
                    end
                    band_mean = mean(band_tc(task_mask), 'omitnan');
                    matchedBlocks(r) = matchedBlocks(r) + 1;
                else
                    idxMissing = t >= onset + analysis_window_sec(1) & t <= onset + analysis_window_sec(2);
                    missingBlockMask(idxMissing) = true;
                end
                alignRows(end+1,:) = {subName, band_name, roi_names{r}, blk, onset, dur, eegAvailable, eegCount, band_mean}; %#ok<AGROW>
            end
        end

        xEegBand = zeros(nT, nRoi);
        for r = 1:nRoi
            x = conv(uEegBand(:,r), hrf);
            xEegBand(:,r) = x(1:nT);
        end

        if zscoreEegRegressors
            for r = 1:nRoi
                reg = xEegBand(:,r);
                mask = isfinite(reg) & reg ~= 0;
                if any(mask)
                    reg_mean = mean(reg(mask), 'omitnan');
                    reg_sd = std(reg(mask), 'omitnan');
                    if ~isfinite(reg_sd) || reg_sd <= 0
                        reg_sd = 1;
                    end
                    xEegBand(:,r) = (reg - reg_mean) / reg_sd;
                end
            end
        end

        if orthogonalizeEegToTask && includeTaskRegressors
            Xtask = xTask;
            if rank(Xtask) > 0
                for r = 1:nRoi
                    xEegBand(:,r) = xEegBand(:,r) - Xtask * (Xtask \ xEegBand(:,r));
                end
            end
        end

        %% Design matrix
        X = [];
        regNames = {};
        if includeTaskRegressors
            X = [X xTask];
            regNames = [regNames {'task'}];
        end
        if includeEegRegressors
            for r = 1:nRoi
                X = [X xEegBand(:,r)];
                regNames = [regNames {sprintf('eeg_%s_%s', band_tag, roi_names{r})}];
            end
        end
        if isempty(X)
            fprintf('[%s] WARNING: No regressors selected for %s (%s). Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName, band_name);
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
            if useRoiSpecificOnly
                roiFromTerm = '';
                if strncmp(termName, 'eeg_', 4)
                    parts = strsplit(termName, '_');
                    if numel(parts) >= 3
                        roiFromTerm = strjoin(parts(3:end), '_');
                    end
                end
                if ~isempty(roiFromTerm)
                    if ~strcmp(roiLabel(mi), roiFromTerm)
                        continue;
                    end
                end
            end
            rowCell(end+1,:) = {subName, taskName, mi, src(mi), det(mi), char(chrom(mi)), char(roiLabel(mi)), ...
                termName, beta(ri,mi), seReg(ri,mi), tReg(ri,mi), pReg(ri,mi)}; %#ok<AGROW>
        end
    end

    T = cell2table(rowCell, 'VariableNames', ...
        {'subjId','task','measIndex','src','det','chromophore','roi','term','beta','SE','t','p'});

    outMat = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm_%s.mat', subName, taskName, band_tag));
    outCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm_table_%s.csv', subName, taskName, band_tag));
    save(outMat, 'beta', 'seReg', 'tReg', 'pReg', 'dof', 'regNames', 'X', ...
        'chrom', 'src', 'det', 'roiLabel', 'roiMask', 'glmMask', 'qcMaskAll', ...
        'eeg_times_ref', 'analysis_ms_use', 'task_window_ms_use', 'band_defs', ...
        'band_name', 'band_range', 'band_count', 'missingBlockMask', '-v7.3');
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
                    if useRoiSpecificOnly
                        roiFromTerm = '';
                        if strncmp(termName, 'eeg_', 4)
                            parts = strsplit(termName, '_');
                            if numel(parts) >= 3
                                roiFromTerm = strjoin(parts(3:end), '_');
                            end
                        end
                        if ~isempty(roiFromTerm)
                            if ~strcmp(roiName, roiFromTerm)
                                continue;
                            end
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
        outRoiCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_glm_roi_table_%s.csv', subName, taskName, band_tag));
        writetable(roiSummaryTable, outRoiCsv);
        save(outMat, 'roiSummaryTable', '-append');
    end

    %% Alignment table
    alignTable = cell2table(alignRows, 'VariableNames', ...
        {'Subject','Band','ROI','BlockIndex','FnirsOnsetSec','FnirsDurSec','EegBlockAvailable','EegGoodTrials','EegBandMean'});
    alignCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_fnirs_alignment_%s.csv', subName, taskName, band_tag));
    writetable(alignTable, alignCsv);

    %% Figures
    if saveFigures
        nPanels = 0;
        if includeTaskRegressors
            nPanels = nPanels + 1;
        end
        if includeEegRegressors
            nPanels = nPanels + 1;
        end
        if nPanels > 0
            figH = figure('Visible','off');
            panelIdx = 0;
            if includeTaskRegressors
                panelIdx = panelIdx + 1;
                subplot(nPanels,1,panelIdx);
                plot(t, xTask, 'Color', [0.2 0.6 0.9]);
                title('Task regressor', 'Interpreter','none');
                legend({'task'}, 'Location','best');
            end
            if includeEegRegressors
                panelIdx = panelIdx + 1;
                subplot(nPanels,1,panelIdx);
                hold on;
                for r = 1:nRoi
                    plot(t, xEegBand(:,r));
                end
                hold off;
                title(sprintf('EEG %s regressors (HRF-convolved)', band_name), 'Interpreter','none');
                legend(roi_names, 'Location','best');
            end
            xlabel('Time (s)');
            figPath = fullfile(subOutDir, sprintf('%s_task-%s_eeg_informed_regressors_%s.png', subName, taskName, band_tag));
            saveas(figH, figPath);
            close(figH);
        end
    end

    %% Summary rows
    missingPct = 0;
    if excludeMissingEegBlocks && includeEegRegressors
        missingPct = 100 * mean(missingBlockMask);
    end
    for r = 1:nRoi
        summaryRows(end+1,:) = {subName, band_name, roi_names{r}, nBlocks_eeg, blocks_good(r,bi), nTaskFnirs, matchedBlocks(r), missingPct}; %#ok<AGROW>
    end
    end

    if processOnlyFirst
        break;
    end
end

%% Save summary
if ~isempty(summaryRows)
    summaryTable = cell2table(summaryRows, 'VariableNames', ...
        {'Subject','Band','ROI','EEG_Blocks','EEG_GoodBlocks','FNIRS_Blocks','MatchedBlocks','MissingTimePct'});
    summaryCsv = fullfile(outRoot, 'STEP2_eeg_informed_fnirs_summary.csv');
    writetable(summaryTable, summaryCsv);
end

fprintf('[%s] STEP2 EEG-informed fNIRS complete. Outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), outRoot);
