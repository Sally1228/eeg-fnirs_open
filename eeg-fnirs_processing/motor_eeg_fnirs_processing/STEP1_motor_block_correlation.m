% EEG-fNIRS block-level correlation for motor task (left/right hand).
% Assumptions & Notes:
% - EEG epochs are time-locked to block start markers (1=left, 2=right).
% - EEG data file naming follows *_task-motor_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set
%   (fallback patterns are checked below).
% - fNIRS block averages are from STEP3_BlockAvarage and saved to *_task-motor_blockavg_table.csv.
% - Block order is derived from EEG event latency; fNIRS blockIndex is chronological within condition.
% - If blocks are missing in either modality, use the intersection of block indices.
% - Correlations are Pearson; p-values require Statistics Toolbox (tcdf).

%% Config
close all; clearvars; clc;

taskName = 'motor';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset patterns (epoched)
eegSetPatterns = {'%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'};

% fNIRS block average table
fnirsBlockTablePattern = '%s_task-%s_blockavg_table.csv';

% EEG condition mapping
cond_labels = {'left', 'right'};
cond_markers = [1 2];

% EEG ROI labels
roi_left_labels = {'C3','C1','C5','CP3','FC3'};
roi_right_labels = {'C4','C2','C6','FC4','CP4'};

% fNIRS ROI names (contralateral mapping)
fnirs_roi_left = 'LeftCentral';
fnirs_roi_right = 'RightCentral';

% EEG time windows
baseline_ms = [-2000 0];
task_window_ms = [0 30000];

% EEG bands
mu_band = [8 12];
beta_band = [13 30];
min_samples = 32;

% Correlation settings
minBlocksForCorr = 3;

%% Paths
% script_dir = '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg_fnirs_processing/motor_eeg_fnirs_processing'
script_dir = fileparts(mfilename('fullpath'));

repo_dir = fileparts(fileparts(script_dir));
eeg_root = fullfile(repo_dir, 'eeg', 'motor');
fnirs_root = fullfile(repo_dir, 'fnirs', 'motor');
out_root = fullfile(repo_dir, 'eeg_fnirs_processing', 'motor_eeg_fnirs_processing','output', 'STEP1_motor_block_correlation');
if exist(out_root, 'dir') ~= 7
    mkdir(out_root);
end

fprintf('[%s] Assumptions: markers {1=left,2=right}, block order by EEG event latency, fNIRS blockIndex from STEP3.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% Discover subjects
if isempty(subjList)
    eegDirs = dir(fullfile(eeg_root, 'sub-*'));
    eegList = {eegDirs([eegDirs.isdir]).name};
    fnirsDirs = dir(fullfile(fnirs_root, 'sub-*'));
    fnirsList = {fnirsDirs([fnirsDirs.isdir]).name};
    subjList = intersect(eegList, fnirsList, 'stable');
end
if isempty(subjList)
    error('No common subjects found between EEG and fNIRS folders.');
end

fprintf('[%s] Subjects: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), strjoin(subjList, ', '));

%% Init containers
blockRows = cell(0, 7);
alignRows = cell(0, 6);

hasTcdf = exist('tcdf', 'file') == 2;
if ~hasTcdf
    fprintf('[%s] WARNING: tcdf not found; p-values will be NaN.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

%% Main loop
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>

for si = 1:numel(subjList)
    subName = subjList{si};
    fprintf('\n[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        subName, si, numel(subjList));

    % --- locate EEG file ---
    eegFile = '';
    for p = 1:numel(eegSetPatterns)
        candidate = fullfile(eeg_root, subName, sprintf(eegSetPatterns{p}, subName, taskName));
        if exist(candidate, 'file') == 2
            eegFile = candidate;
            break;
        end
    end
    if isempty(eegFile)
        fprintf('[%s] WARNING: EEG .set not found for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % --- locate fNIRS block table ---
    fnirsFile = fullfile(fnirs_root, subName, 'blockavg', ...
        sprintf(fnirsBlockTablePattern, subName, taskName));
    if exist(fnirsFile, 'file') ~= 2
        tmp = dir(fullfile(fnirs_root, subName, 'blockavg', '*_blockavg_table.csv'));
        if isempty(tmp)
            tmp = dir(fullfile(fnirs_root, subName, '*_blockavg_table.csv'));
        end
        if ~isempty(tmp)
            fnirsFile = fullfile(tmp(1).folder, tmp(1).name);
        else
            fnirsFile = '';
        end
    end
    if isempty(fnirsFile)
        fprintf('[%s] WARNING: fNIRS block table not found for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % --- load data ---
    [eegPath, eegName, eegExt] = fileparts(eegFile);
    EEG = pop_loadset('filename', [eegName eegExt], 'filepath', eegPath);

    if ~isfield(EEG, 'epoch') || EEG.trials == 0
        fprintf('[%s] WARNING: no epochs for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    fnirsT = readtable(fnirsFile);
    reqCols = {'condition','roi','chromophore','blockIndex','blockMean'};
    if ~all(ismember(reqCols, fnirsT.Properties.VariableNames))
        fprintf('[%s] WARNING: fNIRS table missing required columns for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % --- time windows ---
    epoch_start_ms = EEG.xmin * 1000;
    epoch_end_ms = EEG.xmax * 1000;
    baseline_ms_use = [max(baseline_ms(1), epoch_start_ms), min(baseline_ms(2), min(0, epoch_end_ms))];
    task_ms_use = [max(task_window_ms(1), epoch_start_ms), min(task_window_ms(2), epoch_end_ms)];
    if baseline_ms_use(2) <= baseline_ms_use(1) || task_ms_use(2) <= task_ms_use(1)
        fprintf('[%s] WARNING: invalid baseline/task window for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    times_ms = [];
    if isfield(EEG, 'times') && ~isempty(EEG.times)
        times_ms = EEG.times;
    else
        times_ms = linspace(epoch_start_ms, epoch_end_ms, EEG.pnts);
    end
    baseline_mask = times_ms >= baseline_ms_use(1) & times_ms <= baseline_ms_use(2);
    task_mask = times_ms >= task_ms_use(1) & times_ms <= task_ms_use(2);

    if sum(baseline_mask) < min_samples || sum(task_mask) < min_samples
        fprintf('[%s] WARNING: not enough samples for baseline/task in %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % --- ROI indices ---
    chan_labels = {EEG.chanlocs.labels};
    chan_labels = cellfun(@(x) upper(strtrim(x)), chan_labels, 'UniformOutput', false);
    roi_left_upper = cellfun(@upper, roi_left_labels, 'UniformOutput', false);
    roi_right_upper = cellfun(@upper, roi_right_labels, 'UniformOutput', false);
    left_idx = find(ismember(chan_labels, roi_left_upper));
    right_idx = find(ismember(chan_labels, roi_right_upper));
    if isempty(left_idx) || isempty(right_idx)
        fprintf('[%s] WARNING: EEG ROI channels missing for %s. Skipping subject.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    % --- good trials mask ---
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

    % --- determine bin and event latency per epoch ---
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

    %% Per-condition processing
    for c = 1:numel(cond_labels)
        condName = cond_labels{c};
        condMarker = cond_markers(c);

        if c == 1
            roi_idx = right_idx; % left hand -> right ROI
            fnirs_roi = fnirs_roi_right;
        else
            roi_idx = left_idx; % right hand -> left ROI
            fnirs_roi = fnirs_roi_left;
        end

        cond_all_trials = find(epoch_bin == condMarker);
        cond_good_trials = find(good_trials & epoch_bin == condMarker);

        % Block index mapping by latency (all trials)
        cond_block_index = nan(1, EEG.trials);
        if ~isempty(cond_all_trials)
            cond_lat = epoch_event_latency(cond_all_trials);
            valid_lat = ~isnan(cond_lat);
            sorted_trials = [];
            if any(valid_lat)
                [~, sort_idx] = sort(cond_lat(valid_lat));
                sorted_trials = cond_all_trials(valid_lat);
                sorted_trials = sorted_trials(sort_idx);
            end
            if any(~valid_lat)
                sorted_trials = [sorted_trials, cond_all_trials(~valid_lat)];
            end
            cond_block_index(sorted_trials) = 1:numel(sorted_trials);
        end

        % EEG block features for good trials
        eeg_block_idx = [];
        eeg_mu_erd = [];
        eeg_beta_erd = [];

        for ti = 1:numel(cond_good_trials)
            tr = cond_good_trials(ti);
            block_idx = cond_block_index(tr);
            if isnan(block_idx)
                continue;
            end

            roi_data = squeeze(mean(EEG.data(roi_idx, :, tr), 1));
            roi_data = double(roi_data(:));

            base_data = roi_data(baseline_mask);
            task_data = roi_data(task_mask);
            base_data = base_data(isfinite(base_data));
            task_data = task_data(isfinite(task_data));

            if numel(base_data) < min_samples || numel(task_data) < min_samples
                continue;
            end

            base_data = base_data - mean(base_data);
            task_data = task_data - mean(task_data);

            % Baseline PSD
            Nbase = numel(base_data);
            nfft_base = 2^nextpow2(Nbase);
            Xb = fft(base_data, nfft_base);
            Pxx_base = (abs(Xb).^2) / (EEG.srate * Nbase);
            Pxx_base = Pxx_base(1:floor(nfft_base/2)+1);
            if numel(Pxx_base) > 2
                Pxx_base(2:end-1) = 2 * Pxx_base(2:end-1);
            end
            f_base = (0:floor(nfft_base/2)) * (EEG.srate / nfft_base);
            mu_idx_base = f_base >= mu_band(1) & f_base <= mu_band(2);
            beta_idx_base = f_base >= beta_band(1) & f_base <= beta_band(2);
            mu_base = trapz(f_base(mu_idx_base), Pxx_base(mu_idx_base));
            beta_base = trapz(f_base(beta_idx_base), Pxx_base(beta_idx_base));

            % Task PSD
            Ntask = numel(task_data);
            nfft_task = 2^nextpow2(Ntask);
            Xt = fft(task_data, nfft_task);
            Pxx_task = (abs(Xt).^2) / (EEG.srate * Ntask);
            Pxx_task = Pxx_task(1:floor(nfft_task/2)+1);
            if numel(Pxx_task) > 2
                Pxx_task(2:end-1) = 2 * Pxx_task(2:end-1);
            end
            f_task = (0:floor(nfft_task/2)) * (EEG.srate / nfft_task);
            mu_idx_task = f_task >= mu_band(1) & f_task <= mu_band(2);
            beta_idx_task = f_task >= beta_band(1) & f_task <= beta_band(2);
            mu_task = trapz(f_task(mu_idx_task), Pxx_task(mu_idx_task));
            beta_task = trapz(f_task(beta_idx_task), Pxx_task(beta_idx_task));

            mu_erd = NaN;
            beta_erd = NaN;
            if isfinite(mu_base) && mu_base > 0 && isfinite(mu_task)
                mu_erd = (mu_task / mu_base - 1) * 100;
            end
            if isfinite(beta_base) && beta_base > 0 && isfinite(beta_task)
                beta_erd = (beta_task / beta_base - 1) * 100;
            end

            eeg_block_idx(end+1, 1) = block_idx; %#ok<AGROW>
            eeg_mu_erd(end+1, 1) = mu_erd; %#ok<AGROW>
            eeg_beta_erd(end+1, 1) = beta_erd; %#ok<AGROW>
        end

        % fNIRS block means (contralateral ROI)
        condMask = strcmpi(fnirsT.condition, condName) & strcmpi(fnirsT.roi, fnirs_roi);
        hboMask = condMask & strcmpi(fnirsT.chromophore, 'HbO');
        hbrMask = condMask & strcmpi(fnirsT.chromophore, 'HbR');

        block_idx_hbo = fnirsT.blockIndex(hboMask);
        block_idx_hbr = fnirsT.blockIndex(hbrMask);
        hbo_vals = fnirsT.blockMean(hboMask);
        hbr_vals = fnirsT.blockMean(hbrMask);

        if ~isempty(block_idx_hbo)
            [hbo_unique, ~, hbo_grp] = unique(block_idx_hbo);
            hbo_vals = accumarray(hbo_grp, hbo_vals, [], @(x) mean(x, 'omitnan'));
            block_idx_hbo = hbo_unique;
        end
        if ~isempty(block_idx_hbr)
            [hbr_unique, ~, hbr_grp] = unique(block_idx_hbr);
            hbr_vals = accumarray(hbr_grp, hbr_vals, [], @(x) mean(x, 'omitnan'));
            block_idx_hbr = hbr_unique;
        end

        [fnirs_block_idx, idx_hbo, idx_hbr] = intersect(block_idx_hbo, block_idx_hbr);
        fnirs_hbo = hbo_vals(idx_hbo);
        fnirs_hbr = hbr_vals(idx_hbr);

        % Align EEG and fNIRS blocks
        [common_blocks, idx_fnirs, idx_eeg] = intersect(fnirs_block_idx, eeg_block_idx);

        alignRows(end+1, :) = {subName, condName, numel(cond_all_trials), numel(eeg_block_idx), ...
            numel(fnirs_block_idx), numel(common_blocks)}; %#ok<AGROW>

        for bi = 1:numel(common_blocks)
            blockRows(end+1, :) = {subName, condName, common_blocks(bi), ...
                eeg_mu_erd(idx_eeg(bi)), eeg_beta_erd(idx_eeg(bi)), ...
                fnirs_hbo(idx_fnirs(bi)), fnirs_hbr(idx_fnirs(bi))}; %#ok<AGROW>
        end
    end
end

%% Build tables
if isempty(blockRows)
    error('No matched blocks found. Check EEG/fNIRS inputs and ROI/marker settings.');
end

blockTable = cell2table(blockRows, 'VariableNames', ...
    {'Subject','Condition','BlockIndex','EEG_MuERD_Pct','EEG_BetaERD_Pct','FNIRS_HbO_BlockMean','FNIRS_HbR_BlockMean'});

alignmentTable = cell2table(alignRows, 'VariableNames', ...
    {'Subject','Condition','EEG_TotalBlocks','EEG_GoodBlocks','FNIRS_Blocks','MatchedBlocks'});

%% Correlation analysis
subjects = unique(blockTable.Subject, 'stable');
analysisNames = {'Left', 'Right', 'All'};
bandNames = {'Mu', 'Beta'};
eegCols = {'EEG_MuERD_Pct', 'EEG_BetaERD_Pct'};
chromNames = {'HbO', 'HbR'};
fnirsCols = {'FNIRS_HbO_BlockMean', 'FNIRS_HbR_BlockMean'};

corrRows = cell(0, 8);
for si = 1:numel(subjects)
    subName = subjects{si};
    for ai = 1:numel(analysisNames)
        analysisName = analysisNames{ai};
        if strcmpi(analysisName, 'Left')
            condMask = strcmpi(blockTable.Condition, 'left');
        elseif strcmpi(analysisName, 'Right')
            condMask = strcmpi(blockTable.Condition, 'right');
        else
            condMask = true(height(blockTable), 1);
        end
        subMask = strcmp(blockTable.Subject, subName);
        mask = subMask & condMask;
        if ~any(mask)
            continue;
        end

        for bi = 1:numel(bandNames)
            x = blockTable.(eegCols{bi})(mask);
            for ci = 1:numel(chromNames)
                y = blockTable.(fnirsCols{ci})(mask);
                valid = isfinite(x) & isfinite(y);
                n = sum(valid);
                r = NaN;
                p = NaN;
                detail = '';
                if n >= minBlocksForCorr && std(x(valid)) > 0 && std(y(valid)) > 0
                    R = corrcoef(x(valid), y(valid));
                    r = R(1, 2);
                    if hasTcdf
                        tVal = r * sqrt((n - 2) / max(eps, (1 - r^2)));
                        p = 2 * (1 - tcdf(abs(tVal), n - 2));
                    end
                    block_ids = blockTable.BlockIndex(mask);
                    block_ids = block_ids(valid);
                    detail = strjoin(arrayfun(@(v) sprintf('%d', v), block_ids(:)', 'UniformOutput', false), ';');
                end
                corrRows(end+1, :) = {subName, analysisName, bandNames{bi}, chromNames{ci}, n, r, p, detail}; %#ok<AGROW>
            end
        end
    end
end

% Group-level correlation across all subjects
for ai = 1:numel(analysisNames)
    analysisName = analysisNames{ai};
    if strcmpi(analysisName, 'Left')
        condMask = strcmpi(blockTable.Condition, 'left');
    elseif strcmpi(analysisName, 'Right')
        condMask = strcmpi(blockTable.Condition, 'right');
    else
        condMask = true(height(blockTable), 1);
    end
    idx_all = find(condMask);
    for bi = 1:numel(bandNames)
        x = blockTable.(eegCols{bi})(condMask);
        for ci = 1:numel(chromNames)
            y = blockTable.(fnirsCols{ci})(condMask);
            valid = isfinite(x) & isfinite(y);
            n = sum(valid);
            r = NaN;
            p = NaN;
            detail = '';
            if n >= minBlocksForCorr && std(x(valid)) > 0 && std(y(valid)) > 0
                R = corrcoef(x(valid), y(valid));
                r = R(1, 2);
                if hasTcdf
                    tVal = r * sqrt((n - 2) / max(eps, (1 - r^2)));
                    p = 2 * (1 - tcdf(abs(tVal), n - 2));
                end
                idx_valid = idx_all(valid);
                detail = strjoin(unique(blockTable.Subject(idx_valid), 'stable'), ';');
            end
            corrRows(end+1, :) = {'GroupAll', analysisName, bandNames{bi}, chromNames{ci}, n, r, p, detail}; %#ok<AGROW>
        end
    end
end

corrTable = cell2table(corrRows, 'VariableNames', ...
    {'Subject','Analysis','EEGBand','FNIRSChrom','N','R','P','Detail'});

%% Save outputs
blockCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_features.csv', taskName));
alignCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_alignment.csv', taskName));
corrCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_correlation.csv', taskName));

writetable(blockTable, blockCsv);
writetable(alignmentTable, alignCsv);
writetable(corrTable, corrCsv);

save(fullfile(out_root, sprintf('%s_eeg_fnirs_block_results.mat', taskName)), ...
    'blockTable', 'alignmentTable', 'corrTable', 'baseline_ms', 'task_window_ms', ...
    'mu_band', 'beta_band', 'roi_left_labels', 'roi_right_labels');

fprintf('[%s] Done. Outputs in %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), out_root);
