% STEP7 Block-level band power and ERSP analysis for cate task
% Operates on epoched data from STEP6 (artifact marked, not removed)
% Outputs Excel tables and figures in outputs/STEP7_power

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Each epoch is one 60 s block time-locked to marker=1.
% - Epoch range is about -2000 to 60000 ms (baseline + task).
% - Baseline window is [-2000 0] ms; task window is [0 60000] ms.
% - Rest periods between blocks are not analyzed in this script.
% - Good blocks are those not marked in EEG.reject and without boundary events.
% - ROI labels follow standard 10-5 naming; missing channels are skipped.

fprintf('[%s] Assumption: each epoch is one task block (marker=1).\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Assumption: baseline [-2000 0] ms, task [0 60000] ms.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Assumption: rest periods between blocks are not analyzed here.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
data_dir = fileparts(project_dir);
output_dir = fullfile(project_dir, 'outputs', 'STEP7_power');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

SUB = {'sub-001','sub-002','sub-003','sub-004','sub-005','sub-008','sub-009'};
task = 'cate';

block_labels = {'Vegetable','Animal','Flower','Instrument','Fruit','Appliance'};
nBlocks_expected = 6;

roi_names = {'LeftFrontal','LeftTemporal','RightFrontal','RightTemporal'};
roi_labels = {
    {'F1','F3','F5','AF3','FC1','FC3','FC5','C1','C3'}, ...
    {'FT7','T7','TP7','CP5','C5'}, ...
    {'F2','F4','F6','AF4','FC2','FC4','FC6','C2','C4'}, ...
    {'FT8','T8','TP8','CP6','C6'} ...
};

band_defs = {
    'Delta', [1 4];
    'Theta', [4 8];
    'Alpha', [8 13];
    'Beta',  [13 30]
};

baseline_ms = [-2000 0];
task_window_ms = [0 60000];
analysis_ms = [-2000 60000];

freq_range = [1 40];
tfr_freqs = 1:1:40;
tfr_cycles = [3 0.5];
tfr_timesout = 200;
tfr_padratio = 2;

pwelch_window_sec = 2;
pwelch_overlap_pct = 0.5;

fprintf('[%s] STEP7 power/ERSP analysis start\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% Initialize
nSub = numel(SUB);
nRoi = numel(roi_names);
nBand = size(band_defs, 1);

bandpower_all = nan(nSub, nBlocks_expected, nRoi, nBand);
bandpower_rows = {};
ersp_summary_rows = {};

times_ref = [];
freqs_ref = [];
ersp_all = [];

use_pwelch = exist('pwelch', 'file') == 2;
use_spectopo = exist('spectopo', 'file') == 2;

if ~use_pwelch && ~use_spectopo
    fprintf('[%s] Warning: neither pwelch nor spectopo found. Band power may be skipped.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

if exist('newtimef', 'file') ~= 2
    fprintf('[%s] Warning: newtimef not found. ERSP will be skipped.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

%% Main loop
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for s = 1:nSub
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        SUB{s}, s, nSub);

    Subject_Path = fullfile(data_dir, SUB{s});

    set_candidates = {
        [SUB{s} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'], ...
        [SUB{s} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp.set'], ...
        [SUB{s} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch.set']
    };

    set_file = '';
    for k = 1:numel(set_candidates)
        if exist(fullfile(Subject_Path, set_candidates{k}), 'file')
            set_file = set_candidates{k};
            break;
        end
    end

    if isempty(set_file)
        fprintf('[%s] Warning: no epoched set file found for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

    EEG = pop_loadset('filename', set_file, 'filepath', Subject_Path);

    if ~isfield(EEG, 'epoch') || EEG.trials == 0
        fprintf('[%s] Warning: no epochs for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

    epoch_start_ms = EEG.times(1);
    epoch_end_ms = EEG.times(end);

    baseline_ms_use = [max(baseline_ms(1), epoch_start_ms), min(baseline_ms(2), min(0, epoch_end_ms))];
    analysis_ms_use = [max(analysis_ms(1), epoch_start_ms), min(analysis_ms(2), epoch_end_ms)];
    task_window_ms_use = [max(task_window_ms(1), epoch_start_ms), min(task_window_ms(2), epoch_end_ms)];

    if baseline_ms_use(2) <= baseline_ms_use(1)
        fprintf('[%s] Warning: baseline window invalid for %s. Using full pre-zero range.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        baseline_ms_use = [epoch_start_ms, min(0, epoch_end_ms)];
    end

    if analysis_ms_use(2) <= analysis_ms_use(1)
        fprintf('[%s] Warning: analysis window invalid for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

    task_idx = EEG.times >= task_window_ms_use(1) & EEG.times <= task_window_ms_use(2);

    % Build good trial mask from reject flags + boundary events
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

    if sum(good_trials) == 0
        fprintf('[%s] Warning: no good blocks for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

    chan_labels = {EEG.chanlocs.labels};
    chan_labels = cellfun(@(x) upper(strtrim(x)), chan_labels, 'UniformOutput', false);

    roi_idx_all = cell(1, nRoi);
    for r = 1:nRoi
        roi_upper = cellfun(@upper, roi_labels{r}, 'UniformOutput', false);
        roi_idx_all{r} = find(ismember(chan_labels, roi_upper));
        if isempty(roi_idx_all{r})
            fprintf('[%s] Warning: ROI %s channels not found for %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), roi_names{r}, SUB{s});
        end
    end

    nBlocks_use = min(nBlocks_expected, EEG.trials);
    if EEG.trials ~= nBlocks_expected
        fprintf('[%s] Note: %s has %d epochs (expected %d). Using first %d.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s}, EEG.trials, nBlocks_expected, nBlocks_use);
    end

    for b = 1:nBlocks_use
        if ~good_trials(b)
            fprintf('[%s] Block %d for %s marked bad. Skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), b, SUB{s});
            continue;
        end

        for r = 1:nRoi
            roi_idx = roi_idx_all{r};
            if isempty(roi_idx)
                continue;
            end

            roi_epoch = squeeze(mean(EEG.data(roi_idx, :, b), 1));
            if isvector(roi_epoch)
                roi_epoch = roi_epoch(:);
            end

            % Band power on task window
            if any(task_idx) && (use_pwelch || use_spectopo)
                roi_task = roi_epoch(task_idx);
                roi_task = roi_task(:);

                if use_pwelch
                    win = max(2, round(pwelch_window_sec * EEG.srate));
                    noverlap = round(win * pwelch_overlap_pct);
                    nfft = max(256, 2^nextpow2(win));
                    [pxx, f] = pwelch(double(roi_task), win, noverlap, nfft, EEG.srate);
                else
                    [pxx_db, f] = spectopo(double(roi_task') , 0, EEG.srate, ...
                        'freqrange', freq_range, 'plot', 'off');
                    pxx = 10 .^ (pxx_db / 10);
                    f = f(:);
                end

                for bd = 1:nBand
                    band_name = band_defs{bd, 1};
                    band_range = band_defs{bd, 2};
                    band_mask = f >= band_range(1) & f < band_range(2);
                    if any(band_mask)
                        band_power = trapz(f(band_mask), pxx(band_mask));
                    else
                        band_power = NaN;
                    end
                    bandpower_all(s, b, r, bd) = band_power;

                    row_idx = size(bandpower_rows, 1) + 1;
                    bandpower_rows{row_idx, 1} = SUB{s};
                    bandpower_rows{row_idx, 2} = b;
                    if numel(block_labels) == nBlocks_expected
                        bandpower_rows{row_idx, 3} = block_labels{b};
                    else
                        bandpower_rows{row_idx, 3} = '';
                    end
                    bandpower_rows{row_idx, 4} = roi_names{r};
                    bandpower_rows{row_idx, 5} = band_name;
                    bandpower_rows{row_idx, 6} = band_power;
                    bandpower_rows{row_idx, 7} = numel(roi_idx);
                    bandpower_rows{row_idx, 8} = task_window_ms_use(1);
                    bandpower_rows{row_idx, 9} = task_window_ms_use(2);
                end
            end

            % ERSP via newtimef on this block
            if exist('newtimef', 'file') == 2
                [ersp, ~, ~, times, freqs] = newtimef(double(roi_epoch), EEG.pnts, ...
                    [epoch_start_ms epoch_end_ms], EEG.srate, tfr_cycles, ...
                    'baseline', baseline_ms_use, 'freqs', tfr_freqs, ...
                    'timesout', tfr_timesout, 'padratio', tfr_padratio, ...
                    'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');

                if isempty(times_ref)
                    times_ref = times;
                    freqs_ref = freqs;
                    ersp_all = nan(nSub, nBlocks_expected, nRoi, numel(freqs_ref), numel(times_ref));
                end

                if numel(times) ~= numel(times_ref) || numel(freqs) ~= numel(freqs_ref) || ...
                        any(abs(times - times_ref) > 1e-6) || any(abs(freqs - freqs_ref) > 1e-6)
                    [t_grid, f_grid] = meshgrid(times, freqs);
                    [t_ref, f_ref] = meshgrid(times_ref, freqs_ref);
                    ersp = interp2(t_grid, f_grid, ersp, t_ref, f_ref, 'linear', NaN);
                end

                ersp_all(s, b, r, :, :) = ersp;

                time_mask = times_ref >= task_window_ms_use(1) & times_ref <= task_window_ms_use(2);
                for bd = 1:nBand
                    band_name = band_defs{bd, 1};
                    band_range = band_defs{bd, 2};
                    freq_mask = freqs_ref >= band_range(1) & freqs_ref < band_range(2);
                    if any(freq_mask) && any(time_mask)
                        ersp_mean = mean(ersp(freq_mask, time_mask), 'omitnan');
                    else
                        ersp_mean = NaN;
                    end

                    row_idx = size(ersp_summary_rows, 1) + 1;
                    ersp_summary_rows{row_idx, 1} = SUB{s};
                    ersp_summary_rows{row_idx, 2} = b;
                    if numel(block_labels) == nBlocks_expected
                        ersp_summary_rows{row_idx, 3} = block_labels{b};
                    else
                        ersp_summary_rows{row_idx, 3} = '';
                    end
                    ersp_summary_rows{row_idx, 4} = roi_names{r};
                    ersp_summary_rows{row_idx, 5} = band_name;
                    ersp_summary_rows{row_idx, 6} = ersp_mean;
                    ersp_summary_rows{row_idx, 7} = baseline_ms_use(1);
                    ersp_summary_rows{row_idx, 8} = baseline_ms_use(2);
                    ersp_summary_rows{row_idx, 9} = task_window_ms_use(1);
                    ersp_summary_rows{row_idx, 10} = task_window_ms_use(2);
                end
            end
        end
    end
end

%% Save tables
excel_path = fullfile(output_dir, 'STEP7_power_results.xlsx');

if ~isempty(bandpower_rows)
    bandpower_tbl = cell2table(bandpower_rows, 'VariableNames', ...
        {'Subject','Block','BlockLabel','ROI','Band','BandPower','NChannels','TaskStart_ms','TaskEnd_ms'});
    writetable(bandpower_tbl, excel_path, 'Sheet', 'BandPower');
else
    bandpower_tbl = table();
end

if ~isempty(ersp_summary_rows)
    ersp_tbl = cell2table(ersp_summary_rows, 'VariableNames', ...
        {'Subject','Block','BlockLabel','ROI','Band','ERSP_dB','BaselineStart_ms','BaselineEnd_ms','TaskStart_ms','TaskEnd_ms'});
    writetable(ersp_tbl, excel_path, 'Sheet', 'ERSP_Summary');
else
    ersp_tbl = table();
end

save(fullfile(output_dir, 'STEP7_power_results.mat'), ...
    'bandpower_all', 'ersp_all', 'times_ref', 'freqs_ref', ...
    'roi_names', 'roi_labels', 'band_defs', 'baseline_ms', 'task_window_ms');

%% Save per-subject outputs
for s = 1:nSub
    sub_id = SUB{s};
    sub_output_dir = fullfile(data_dir, sub_id, 'outputs', 'STEP7_power');
    if ~exist(sub_output_dir, 'dir')
        mkdir(sub_output_dir);
    end

    sub_excel = fullfile(sub_output_dir, 'STEP7_power_results.xlsx');

    if ~isempty(bandpower_tbl)
        sub_mask = strcmp(string(bandpower_tbl.Subject), sub_id);
        sub_bandpower = bandpower_tbl(sub_mask, :);
        if ~isempty(sub_bandpower)
            writetable(sub_bandpower, sub_excel, 'Sheet', 'BandPower');
        end
    end

    if ~isempty(ersp_tbl)
        sub_mask = strcmp(string(ersp_tbl.Subject), sub_id);
        sub_ersp = ersp_tbl(sub_mask, :);
        if ~isempty(sub_ersp)
            writetable(sub_ersp, sub_excel, 'Sheet', 'ERSP_Summary');
        end
    end

    bandpower_sub = squeeze(bandpower_all(s, :, :, :));
    if ~isempty(ersp_all)
        ersp_sub = squeeze(ersp_all(s, :, :, :, :));
    else
        ersp_sub = [];
    end

    save(fullfile(sub_output_dir, 'STEP7_power_results.mat'), ...
        'bandpower_sub', 'ersp_sub', 'times_ref', 'freqs_ref', ...
        'roi_names', 'roi_labels', 'band_defs', 'baseline_ms', 'task_window_ms');
end

%% Figures: band power across blocks (group mean)
if any(~isnan(bandpower_all(:)))
    bandpower_db = 10 * log10(max(bandpower_all, realmin));
    colors = lines(nBand);

    figure('Color', 'w', 'Position', [100 100 1200 800]);
    for r = 1:nRoi
        subplot(2, 2, r);
        hold on;
        for bd = 1:nBand
            y = squeeze(mean(bandpower_db(:, :, r, bd), 1, 'omitnan'));
            plot(1:nBlocks_expected, y, '-o', 'LineWidth', 1.5, 'Color', colors(bd, :));
        end
        hold off;
        grid on;
        xlim([1 nBlocks_expected]);
        if numel(block_labels) == nBlocks_expected
            set(gca, 'XTick', 1:nBlocks_expected, 'XTickLabel', block_labels);
        else
            set(gca, 'XTick', 1:nBlocks_expected);
        end
        xlabel('Block');
        ylabel('Band Power (dB)');
        title(roi_names{r});
        if r == 1
            legend(band_defs(:,1), 'Location', 'best');
        end
    end
    saveas(gcf, fullfile(output_dir, 'BandPower_Blocks.png'));
    close(gcf);
end

%% Figures: ERSP per block (group mean)
if ~isempty(ersp_all) && ~isempty(times_ref) && ~isempty(freqs_ref)
    ersp_vals = ersp_all(:);
    ersp_vals = ersp_vals(~isnan(ersp_vals));
    if isempty(ersp_vals)
        ersp_clim = [-3 3];
    else
        ersp_clim = prctile(ersp_vals, [5 95]);
        if ersp_clim(1) == ersp_clim(2)
            ersp_clim = ersp_clim + [-1 1];
        end
    end

    for b = 1:nBlocks_expected
        figure('Color', 'w', 'Position', [100 100 1200 800]);
        for r = 1:nRoi
            subplot(2, 2, r);
            ersp_mean = squeeze(mean(ersp_all(:, b, r, :, :), 1, 'omitnan'));
            imagesc(times_ref, freqs_ref, ersp_mean);
            axis xy;
            colormap(jet);
            caxis(ersp_clim);
            colorbar;
            hold on;
            line([0 0], [freqs_ref(1) freqs_ref(end)], 'Color', 'k', 'LineStyle', '--');
            line([task_window_ms(2) task_window_ms(2)], [freqs_ref(1) freqs_ref(end)], 'Color', 'k', 'LineStyle', ':');
            hold off;
            xlabel('Time (ms)');
            ylabel('Freq (Hz)');
            title(roi_names{r});
        end
        if exist('sgtitle', 'file') == 2
            if numel(block_labels) == nBlocks_expected
                sgtitle(['ERSP Block ' num2str(b) ' - ' block_labels{b}]);
            else
                sgtitle(['ERSP Block ' num2str(b)]);
            end
        end
        saveas(gcf, fullfile(output_dir, sprintf('ERSP_Block_%02d.png', b)));
        close(gcf);
    end
end

fprintf('[%s] STEP7 power/ERSP analysis done\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
