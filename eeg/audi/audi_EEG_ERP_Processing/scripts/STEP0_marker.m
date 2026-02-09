% Script #0 Convert markers from trigger channels

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/audi/sub-XXX).

%script_dir = '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/audi/audi_EEG_ERP_Processing/scripts';
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
DIR = fileparts(project_dir);
out_dir = fullfile(project_dir,'outputs');


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), DIR);
fprintf('[%s] Script #0 Convert markers from trigger channels\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-002','sub-003', 'sub-004', 'sub-006','sub-008', 'sub-009'};

% Trigger channel indices in the EEG data 运动任务用 38、39、40通道
trigChans = [38 39 40];

% Options
plot_qc = true;
plot_first_only = true;
keep_boundary_events = true;
export_marker_tables = true; %如果需要输出marker到excel检查，改成true

% Binarization parameters
zero_tol = 1e-6;
kmeans_reps = 5;
min_sep_factor = 5;

% Optional output tables
if export_marker_tables
    all_marker_rows = {'file', 'code', 'time_s'};
    summary_rows = {'file', 'code', 'count'};
end

% Open EEGLAB once
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<NASGU,ASGLU>

task = 'audi';
task_plotted = struct();
task_plotted.(task) = false;

for i = 1:numel(SUB)
    sid = SUB{i};
    subFolder = fullfile(DIR, sid);

    baseName = sprintf('%s_task-audi_eeg', sid);
    inFile = fullfile(subFolder, [baseName '.edf']);
    outFile = [baseName '.set'];
    file_name = [baseName '.edf'];

    fprintf('\n=== (%d/%d) sub-%s ===\n', i, numel(SUB), sid);
    fprintf('Input : %s\n', inFile);

    if ~exist(inFile, 'file')
        warning('EDF not found, skipping: %s', inFile);
        continue;
    end

    EEG = pop_biosig(inFile);

    if max(trigChans) > EEG.nbchan
        warning('Trigger channel index exceeds channel count (nbchan=%d), skipping: %s', ...
            EEG.nbchan, inFile);
        continue;
    end

    trig_labels = arrayfun(@(c) sprintf('ch%d', c), trigChans, 'UniformOutput', false);
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= max(trigChans)
        for ii = 1:numel(trigChans)
            lbl = EEG.chanlocs(trigChans(ii)).labels;
            if ~isempty(lbl)
                trig_labels{ii} = lbl;
            end
        end
    end

    %% Binarize trigger channels
    X = double(EEG.data(trigChans, :)); % [nCh x N]
    fs = EEG.srate;
    N = size(X, 2);
    B = zeros(size(X));

    for ich = 1:size(X, 1)
        x_raw = X(ich, :);
        x_raw(isnan(x_raw)) = 0;

        zero_mask = abs(x_raw) <= zero_tol;
        x_nz = x_raw(~zero_mask);

        chan_label = trig_labels{ich};

        if isempty(x_nz)
            B(ich, :) = 0;
            fprintf('Channel %d (%s): all zeros, no events.\n', trigChans(ich), chan_label);
            continue;
        end

        x_nz_range = max(x_nz) - min(x_nz);
        if x_nz_range <= zero_tol
            b = zeros(size(x_raw));
            if any(zero_mask)
                b(~zero_mask) = 1;
                fprintf('Channel %d (%s): single nonzero level, treating nonzero as high.\n', ...
                    trigChans(ich), chan_label);
            else
                fprintf('Channel %d (%s): constant level, no events.\n', ...
                    trigChans(ich), chan_label);
            end
            B(ich, :) = b;
            continue;
        end

        x = x_nz - median(x_nz);
        b = zeros(size(x_raw));
        noise = mad(x, 1);
        sep = NaN;

        use_kmeans = exist('kmeans', 'file') == 2;
        if use_kmeans
            try
                [idx, C] = kmeans(x', 2, 'Replicates', kmeans_reps);
                [~, low_idx] = min(C);
                [~, high_idx] = max(C);
                sep = C(high_idx) - C(low_idx);

                min_sep = max(min_sep_factor * noise, zero_tol);
                if sep >= min_sep
                    nz_idx = find(~zero_mask);
                    b(nz_idx(idx == high_idx)) = 1;
                else
                    fprintf('Channel %d (%s): no distinct high level (sep=%.3e, noise=%.3e)\n', ...
                        trigChans(ich), chan_label, sep, noise);
                end
            catch
                use_kmeans = false;
            end
        end

        if ~use_kmeans
            low_level = min(x_nz);
            high_level = max(x_nz);
            sep = high_level - low_level;
            min_sep = max(min_sep_factor * noise, zero_tol);

            if sep >= min_sep
                thr = (low_level + high_level) / 2;
                b(~zero_mask) = x_nz > thr;
            else
                fprintf('Channel %d (%s): no distinct high level (sep=%.3e, noise=%.3e)\n', ...
                    trigChans(ich), chan_label, sep, noise);
            end
        end

        B(ich, :) = b;
    end

    %% Combine trigger bits into event code
    weights = 2.^(0:size(B, 1) - 1);
    code = weights * B;

    event_codes = unique(code);
    fprintf('Event codes: ');
    disp(event_codes);

    %% Find event onsets (rising edges from 0)
    code_prev = [0, code(1:end-1)];
    onset_idx = find(code_prev == 0 & code > 0);
    onset_code = code(onset_idx);
    onset_time = (onset_idx - 1) / fs;

    if isempty(onset_idx)
        warning('No event onsets detected in %s.', file_name);
    end

    % Optional export tables
    if export_marker_tables && ~isempty(onset_idx)
        for k = 1:numel(onset_idx)
            all_marker_rows(end+1, :) = {file_name, double(onset_code(k)), double(onset_time(k))}; %#ok<AGROW>
        end

        [u_marker, ~, ic] = unique(onset_code);
        counts = accumarray(ic, 1);
        for m = 1:numel(u_marker)
            summary_rows(end+1, :) = {file_name, double(u_marker(m)), double(counts(m))}; %#ok<AGROW>
        end
    end

    % QC plot
    if plot_qc && (~plot_first_only || ~task_plotted.(task))
        t = (0:N-1) / fs;
        figure('Name', sprintf('%s | %s', file_name, task), 'NumberTitle', 'off');
        subplot(2,1,1);
        plot(t, B');
        xlabel('Time (s)');
        ylabel('bit value');
        legend(trig_labels, 'Location', 'best');
        title('Binarized trigger signals');
        ylim([-0.2 1.2]);

        subplot(2,1,2);
        plot(t, code);
        hold on;
        if ~isempty(onset_idx)
            stem(onset_time, onset_code, 'filled');
        end
        hold off;
        xlabel('Time (s)');
        ylabel('Code');
        title('Trigger code and detected onsets');

        task_plotted.(task) = true;
    end

    %% Write events into EEG.event
    if keep_boundary_events && ~isempty(EEG.event)
        EEG.event = EEG.event(strcmp({EEG.event.type}, 'boundary'));
    else
        EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
    end

    n0 = numel(EEG.event);
    for k = 1:numel(onset_idx)
        EEG.event(n0 + k).type = num2str(onset_code(k));
        EEG.event(n0 + k).latency = onset_idx(k);
        EEG.event(n0 + k).duration = 0;
    end

    EEG = eeg_checkset(EEG, 'eventconsistency');

    %% Save .set
    EEG = pop_saveset(EEG, 'filename', outFile, 'filepath', subFolder);
    fprintf('Saved: %s\n', fullfile(subFolder, outFile));
end

if export_marker_tables
    writecell(all_marker_rows, fullfile(out_dir, 'marker_events.csv'));
    writecell(summary_rows, fullfile(out_dir, 'marker_counts.csv'));
end
