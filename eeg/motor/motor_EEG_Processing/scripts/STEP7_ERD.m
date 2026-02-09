% STEP7 Motor ERD (mu/beta) analysis for block design
% Operates on individual subject data
% Uses the output from Script #6: Artifact_Rejection.m
% This script computes mu/beta ERD time courses for left/right hand blocks,
% compares contralateral vs ipsilateral ROIs, and saves plots + Excel summary.

% 从 STEP6 的已分段 EEG 中计算左右手 block 的 μ/β 频段 ERD，做对侧/同侧比较，并输出时间曲线图 + Excel 汇总。
% 主要流程
% - Config 区：定义被试列表、条件（bin1=左手，bin2=右手）、ROI（C3/C4 周围）、时间窗（baseline/analysis/task）、频段（μ/β）和 newtimef 参数；输出目录为 eeg/motor/motor_EEG_Processing/outputs/ERD.
% - 逐被试加载数据：读取 *_elist_bins_epoch_interp_ar.set，自动把 baseline_ms / analysis_ms / task_window_ms 裁剪到实际 epoch 范围。
% - 好试次筛选：把 EEG.reject 的所有标记 + boundary 事件对应 epoch 排除掉，剩余为 good trials。
% - ROI 映射：按通道标签找左/右 ROI 的通道索引。
% - 确定每个 epoch 的条件：优先读 EEG.EVENTLIST.eventinfo(ev).bini（ERPLAB bin），否则回退到 event type。
% ERD 计算：
% - 先把 ROI 内通道平均成 1 个时序，再用 newtimef 做时频（baseline 归一化为 dB）。
% - 对 μ/β 频段做频率平均，得到时间曲线。
% - 把 dB 转成百分比变化：(10^(dB/10) - 1) * 100。
% - 结果对齐到统一时间轴后保存。
% - 统计汇总：在 task_window_ms 内平均 ERD，分别输出 ROI 平均值 + 对侧/同侧差值与侧化指数（LI）。
% 输出图：
% - μ/β 的对侧 vs 同侧时间曲线（含 95% CI）
% - 对侧-同侧差值时间曲线
% - Excel 汇总：STEP7_ERD_Summary.xlsx
% - 时间序列数据：STEP7_ERD_Timecourses.mat
% - 图像：ERD_Mu_Timecourse.png、ERD_Beta_Timecourse.png、ERD_Mu_ContraMinusIpsi.png、ERD_Beta_ContraMinusIpsi.png（组平均）




%% Config
close all; clearvars;

fprintf('[%s] Assumption: subject folders live in sibling data directory.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Assumption: bins {1:left, 2:right}, baseline [-2000 0] ms.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% script_dir =  '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/motor/motor_EEG_Processing/scripts'
script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);
input_dir = fullfile(project_dir, 'inputs');
output_dir = fullfile(project_dir, 'outputs', 'ERD');
data_dir = fileparts(project_dir);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

SUB = {'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006','sub-008','sub-009'};
task = 'motor';

cond_labels = {'Left', 'Right'};
cond_markers = [1 2];

roi_left_labels = {'C3','C1','C5','CP3','FC3'};
roi_right_labels = {'C4','C2','C6','FC4','CP4'};
roi_names = {'LeftROI', 'RightROI'};

baseline_ms = [-2000 0];
analysis_ms = [-2000 40000]; %分析-2s至40s， 用于截取/绘制ERD 时间曲线的分析区间。不改变epoch
task_window_ms = [0 30000]; %任务窗口0s至30s，计算平均 ERD 值

mu_band = [8 12];
beta_band = [13 30];
tfr_freqs = [4 40]; %时频分析计算的整体频率范围（Hz），决定 newtimef 输出的频率轴。
tfr_cycles = [3 0.5];%小波周期数设置
tfr_timesout = 200;
tfr_padratio = 2;

ci_multiplier = 1.96;

fprintf('[%s] STEP7 ERD analysis start\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% Initialize containers
nSub = length(SUB);
nCond = 2;
nRoi = 2;
times_ref = [];
erd_mu = [];
erd_beta = [];
good_trials_count = zeros(nSub, nCond);

summary_roi_rows = {};
summary_contra_rows = {};

%% Main processing
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for s = 1:nSub
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        SUB{s}, s, nSub);

    Subject_Path = fullfile(data_dir, SUB{s});
    EEG = pop_loadset('filename', [SUB{s} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'], ...
        'filepath', Subject_Path);

    if ~isfield(EEG, 'epoch') || EEG.trials == 0
        fprintf('[%s] Warning: no epochs for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

    % Adjust windows to available epoch range
    epoch_start_ms = EEG.xmin * 1000;
    epoch_end_ms = EEG.xmax * 1000;
    baseline_ms_use = [max(baseline_ms(1), epoch_start_ms), min(baseline_ms(2), min(0, epoch_end_ms))];
    analysis_ms_use = [max(analysis_ms(1), epoch_start_ms), min(analysis_ms(2), epoch_end_ms)];
    task_window_ms_use = [max(task_window_ms(1), epoch_start_ms), min(task_window_ms(2), epoch_end_ms)];

    if baseline_ms_use(2) <= baseline_ms_use(1)
        fprintf('[%s] Warning: baseline window invalid for %s. Using [%d %d] ms.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s}, round(epoch_start_ms), round(min(0, epoch_end_ms)));
        baseline_ms_use = [epoch_start_ms, min(0, epoch_end_ms)];
    end

    if analysis_ms_use(2) <= analysis_ms_use(1)
        fprintf('[%s] Warning: analysis window invalid for %s. Skipping.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
        continue;
    end

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

    % Map ROI channel indices
    chan_labels = {EEG.chanlocs.labels};
    chan_labels = cellfun(@(x) upper(strtrim(x)), chan_labels, 'UniformOutput', false);
    roi_left_upper = cellfun(@upper, roi_left_labels, 'UniformOutput', false);
    roi_right_upper = cellfun(@upper, roi_right_labels, 'UniformOutput', false);
    left_idx = find(ismember(chan_labels, roi_left_upper));
    right_idx = find(ismember(chan_labels, roi_right_upper));

    if isempty(left_idx)
        fprintf('[%s] Warning: left ROI channels not found for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
    end
    if isempty(right_idx)
        fprintf('[%s] Warning: right ROI channels not found for %s.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s});
    end

    % Determine bin per epoch (time-locking event closest to 0 ms)
    epoch_bin = nan(1, EEG.trials);
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
    end

    % Compute ERD for each condition and ROI
    for c = 1:nCond
        cond_trials = find(good_trials & epoch_bin == cond_markers(c));
        good_trials_count(s, c) = numel(cond_trials);
        if isempty(cond_trials)
            fprintf('[%s] Warning: no good trials for %s %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{s}, cond_labels{c});
            continue;
        end

        for r = 1:nRoi
            if r == 1
                roi_idx = left_idx;
            else
                roi_idx = right_idx;
            end
            if isempty(roi_idx)
                continue;
            end

            roi_data = squeeze(mean(EEG.data(roi_idx, :, cond_trials), 1));
            if isvector(roi_data)
                roi_data = roi_data(:);
            end

            [ersp, ~, ~, times, freqs] = newtimef(double(roi_data), EEG.pnts, ...
                [epoch_start_ms epoch_end_ms], EEG.srate, tfr_cycles, ...
                'baseline', baseline_ms_use, 'freqs', tfr_freqs, ...
                'timesout', tfr_timesout, 'padratio', tfr_padratio, ...
                'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');

            time_mask = times >= analysis_ms_use(1) & times <= analysis_ms_use(2);
            mu_idx = freqs >= mu_band(1) & freqs <= mu_band(2);
            beta_idx = freqs >= beta_band(1) & freqs <= beta_band(2);

            mu_db = mean(ersp(mu_idx, time_mask), 1);
            beta_db = mean(ersp(beta_idx, time_mask), 1);
            mu_pct = (10 .^ (mu_db / 10) - 1) * 100;
            beta_pct = (10 .^ (beta_db / 10) - 1) * 100;
            times_sel = times(time_mask);

            if isempty(times_ref)
                times_ref = times_sel;
                nTimes = numel(times_ref);
                erd_mu = nan(nSub, nCond, nRoi, nTimes);
                erd_beta = nan(nSub, nCond, nRoi, nTimes);
            end

            if numel(times_sel) ~= numel(times_ref) || any(abs(times_sel - times_ref) > 1e-6)
                mu_pct = interp1(times_sel, mu_pct, times_ref, 'linear', 'extrap');
                beta_pct = interp1(times_sel, beta_pct, times_ref, 'linear', 'extrap');
            end

            erd_mu(s, c, r, :) = mu_pct;
            erd_beta(s, c, r, :) = beta_pct;
        end
    end
end

if isempty(times_ref)
    fprintf('[%s] No ERD results computed. Check input data and bins.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    return;
end

%% Summary metrics
task_mask = times_ref >= task_window_ms(1) & times_ref <= task_window_ms(2);
if ~any(task_mask)
    fprintf('[%s] Warning: task window outside analysis range. Using full analysis window.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    task_mask = true(size(times_ref));
end
task_window_ms_use = [times_ref(find(task_mask, 1, 'first')), times_ref(find(task_mask, 1, 'last'))];

roi_mean = nan(nSub, nCond, nRoi, 2);
for s = 1:nSub
    for c = 1:nCond
        for r = 1:nRoi
            mu_series = squeeze(erd_mu(s, c, r, :));
            beta_series = squeeze(erd_beta(s, c, r, :));
            roi_mean(s, c, r, 1) = mean(mu_series(task_mask), 'omitnan');
            roi_mean(s, c, r, 2) = mean(beta_series(task_mask), 'omitnan');
        end
    end
end

band_names = {'Mu', 'Beta'};

roi_row = 0;
for s = 1:nSub
    for c = 1:nCond
        for r = 1:nRoi
            for b = 1:2
                roi_row = roi_row + 1;
                summary_roi_rows{roi_row, 1} = SUB{s};
                summary_roi_rows{roi_row, 2} = cond_labels{c};
                summary_roi_rows{roi_row, 3} = roi_names{r};
                summary_roi_rows{roi_row, 4} = band_names{b};
                summary_roi_rows{roi_row, 5} = roi_mean(s, c, r, b);
                summary_roi_rows{roi_row, 6} = good_trials_count(s, c);
                summary_roi_rows{roi_row, 7} = task_window_ms_use(1);
                summary_roi_rows{roi_row, 8} = task_window_ms_use(2);
            end
        end
    end
end

contra_row = 0;
for s = 1:nSub
    for c = 1:nCond
        if c == 1
            contra_roi = 2; % left hand -> right ROI
            ipsi_roi = 1;
        else
            contra_roi = 1; % right hand -> left ROI
            ipsi_roi = 2;
        end
        for b = 1:2
            contra_mean = roi_mean(s, c, contra_roi, b);
            ipsi_mean = roi_mean(s, c, ipsi_roi, b);
            diff_mean = contra_mean - ipsi_mean;
            contra_mag = -contra_mean;
            ipsi_mag = -ipsi_mean;
            if isnan(contra_mag) || isnan(ipsi_mag) || (contra_mag + ipsi_mag) == 0
                li = NaN;
            else
                li = (contra_mag - ipsi_mag) / (contra_mag + ipsi_mag);
            end
            contra_row = contra_row + 1;
            summary_contra_rows{contra_row, 1} = SUB{s};
            summary_contra_rows{contra_row, 2} = cond_labels{c};
            summary_contra_rows{contra_row, 3} = band_names{b};
            summary_contra_rows{contra_row, 4} = contra_mean;
            summary_contra_rows{contra_row, 5} = ipsi_mean;
            summary_contra_rows{contra_row, 6} = diff_mean;
            summary_contra_rows{contra_row, 7} = li;
            summary_contra_rows{contra_row, 8} = good_trials_count(s, c);
            summary_contra_rows{contra_row, 9} = task_window_ms_use(1);
            summary_contra_rows{contra_row, 10} = task_window_ms_use(2);
        end
    end
end

roi_table = cell2table(summary_roi_rows, 'VariableNames', ...
    {'Subject','Condition','ROI','Band','MeanERD_Pct','GoodTrials','TaskStart_ms','TaskEnd_ms'});
contra_table = cell2table(summary_contra_rows, 'VariableNames', ...
    {'Subject','Condition','Band','ContraMeanERD_Pct','IpsiMeanERD_Pct','ContraMinusIpsi_Pct','LateralizationIndex','GoodTrials','TaskStart_ms','TaskEnd_ms'});

summary_file = fullfile(output_dir, 'STEP7_ERD_Summary.xlsx');
writetable(roi_table, summary_file, 'Sheet', 'ROI');
writetable(contra_table, summary_file, 'Sheet', 'ContraIpsi');

save(fullfile(output_dir, 'STEP7_ERD_Timecourses.mat'), ...
    'times_ref', 'erd_mu', 'erd_beta', 'SUB', 'cond_labels', 'roi_names', ...
    'mu_band', 'beta_band', 'baseline_ms', 'analysis_ms', 'task_window_ms');

%% Group plots: mu/beta ERD time courses
time_sec = times_ref / 1000;

% Helper inline for mean/CI
% Contra/Ipsi indices: Left hand -> (R,L), Right hand -> (L,R)

figure_mu = figure('Color', 'w', 'Position', [100 100 1100 450]);
for c = 1:nCond
    if c == 1
        contra_roi = 2;
        ipsi_roi = 1;
    else
        contra_roi = 1;
        ipsi_roi = 2;
    end
    mu_contra = squeeze(erd_mu(:, c, contra_roi, :));
    mu_ipsi = squeeze(erd_mu(:, c, ipsi_roi, :));

    mu_contra_mean = mean(mu_contra, 1, 'omitnan');
    mu_contra_sd = std(mu_contra, 0, 1, 'omitnan');
    mu_contra_n = sum(~isnan(mu_contra), 1);
    mu_contra_sem = mu_contra_sd ./ sqrt(mu_contra_n);
    mu_contra_ci = ci_multiplier * mu_contra_sem;

    mu_ipsi_mean = mean(mu_ipsi, 1, 'omitnan');
    mu_ipsi_sd = std(mu_ipsi, 0, 1, 'omitnan');
    mu_ipsi_n = sum(~isnan(mu_ipsi), 1);
    mu_ipsi_sem = mu_ipsi_sd ./ sqrt(mu_ipsi_n);
    mu_ipsi_ci = ci_multiplier * mu_ipsi_sem;

    subplot(1, 2, c);
    hold on;

    contra_color = [0.2 0.4 0.8];
    ipsi_color = [0.85 0.4 0.2];

    valid_contra = ~isnan(mu_contra_mean);
    if any(valid_contra)
        t = time_sec(valid_contra);
        m = mu_contra_mean(valid_contra);
        ci = mu_contra_ci(valid_contra);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], contra_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', contra_color, 'LineWidth', 2);
    end

    valid_ipsi = ~isnan(mu_ipsi_mean);
    if any(valid_ipsi)
        t = time_sec(valid_ipsi);
        m = mu_ipsi_mean(valid_ipsi);
        ci = mu_ipsi_ci(valid_ipsi);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], ipsi_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', ipsi_color, 'LineWidth', 2);
    end

    xline(0, 'k--');
    yline(0, 'k:');
    xlabel('Time (s)');
    ylabel('Mu ERD (% vs baseline)');
    title([cond_labels{c} ' hand']);
    legend({'Contra CI','Contra','Ipsi CI','Ipsi'}, 'Location', 'best');
    grid on;
    hold off;
end

saveas(figure_mu, fullfile(output_dir, 'ERD_Mu_Timecourse.png'));
close(figure_mu);

figure_beta = figure('Color', 'w', 'Position', [100 100 1100 450]);
for c = 1:nCond
    if c == 1
        contra_roi = 2;
        ipsi_roi = 1;
    else
        contra_roi = 1;
        ipsi_roi = 2;
    end
    beta_contra = squeeze(erd_beta(:, c, contra_roi, :));
    beta_ipsi = squeeze(erd_beta(:, c, ipsi_roi, :));

    beta_contra_mean = mean(beta_contra, 1, 'omitnan');
    beta_contra_sd = std(beta_contra, 0, 1, 'omitnan');
    beta_contra_n = sum(~isnan(beta_contra), 1);
    beta_contra_sem = beta_contra_sd ./ sqrt(beta_contra_n);
    beta_contra_ci = ci_multiplier * beta_contra_sem;

    beta_ipsi_mean = mean(beta_ipsi, 1, 'omitnan');
    beta_ipsi_sd = std(beta_ipsi, 0, 1, 'omitnan');
    beta_ipsi_n = sum(~isnan(beta_ipsi), 1);
    beta_ipsi_sem = beta_ipsi_sd ./ sqrt(beta_ipsi_n);
    beta_ipsi_ci = ci_multiplier * beta_ipsi_sem;

    subplot(1, 2, c);
    hold on;

    contra_color = [0.2 0.6 0.4];
    ipsi_color = [0.8 0.4 0.6];

    valid_contra = ~isnan(beta_contra_mean);
    if any(valid_contra)
        t = time_sec(valid_contra);
        m = beta_contra_mean(valid_contra);
        ci = beta_contra_ci(valid_contra);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], contra_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', contra_color, 'LineWidth', 2);
    end

    valid_ipsi = ~isnan(beta_ipsi_mean);
    if any(valid_ipsi)
        t = time_sec(valid_ipsi);
        m = beta_ipsi_mean(valid_ipsi);
        ci = beta_ipsi_ci(valid_ipsi);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], ipsi_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', ipsi_color, 'LineWidth', 2);
    end

    xline(0, 'k--');
    yline(0, 'k:');
    xlabel('Time (s)');
    ylabel('Beta ERD (% vs baseline)');
    title([cond_labels{c} ' hand']);
    legend({'Contra CI','Contra','Ipsi CI','Ipsi'}, 'Location', 'best');
    grid on;
    hold off;
end

saveas(figure_beta, fullfile(output_dir, 'ERD_Beta_Timecourse.png'));
close(figure_beta);

%% Contra vs ipsi difference time courses
figure_mu_diff = figure('Color', 'w', 'Position', [100 100 1100 450]);
for c = 1:nCond
    if c == 1
        contra_roi = 2;
        ipsi_roi = 1;
    else
        contra_roi = 1;
        ipsi_roi = 2;
    end
    mu_diff = squeeze(erd_mu(:, c, contra_roi, :) - erd_mu(:, c, ipsi_roi, :));
    mu_diff_mean = mean(mu_diff, 1, 'omitnan');
    mu_diff_sd = std(mu_diff, 0, 1, 'omitnan');
    mu_diff_n = sum(~isnan(mu_diff), 1);
    mu_diff_sem = mu_diff_sd ./ sqrt(mu_diff_n);
    mu_diff_ci = ci_multiplier * mu_diff_sem;

    subplot(1, 2, c);
    hold on;
    diff_color = [0.3 0.3 0.3];
    valid_diff = ~isnan(mu_diff_mean);
    if any(valid_diff)
        t = time_sec(valid_diff);
        m = mu_diff_mean(valid_diff);
        ci = mu_diff_ci(valid_diff);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], diff_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', diff_color, 'LineWidth', 2);
    end
    xline(0, 'k--');
    yline(0, 'k:');
    xlabel('Time (s)');
    ylabel('Mu Contra - Ipsi ERD (%)');
    title([cond_labels{c} ' hand']);
    grid on;
    hold off;
end

saveas(figure_mu_diff, fullfile(output_dir, 'ERD_Mu_ContraMinusIpsi.png'));
close(figure_mu_diff);

figure_beta_diff = figure('Color', 'w', 'Position', [100 100 1100 450]);
for c = 1:nCond
    if c == 1
        contra_roi = 2;
        ipsi_roi = 1;
    else
        contra_roi = 1;
        ipsi_roi = 2;
    end
    beta_diff = squeeze(erd_beta(:, c, contra_roi, :) - erd_beta(:, c, ipsi_roi, :));
    beta_diff_mean = mean(beta_diff, 1, 'omitnan');
    beta_diff_sd = std(beta_diff, 0, 1, 'omitnan');
    beta_diff_n = sum(~isnan(beta_diff), 1);
    beta_diff_sem = beta_diff_sd ./ sqrt(beta_diff_n);
    beta_diff_ci = ci_multiplier * beta_diff_sem;

    subplot(1, 2, c);
    hold on;
    diff_color = [0.3 0.3 0.3];
    valid_diff = ~isnan(beta_diff_mean);
    if any(valid_diff)
        t = time_sec(valid_diff);
        m = beta_diff_mean(valid_diff);
        ci = beta_diff_ci(valid_diff);
        fill([t fliplr(t)], [m - ci fliplr(m + ci)], diff_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, m, 'Color', diff_color, 'LineWidth', 2);
    end
    xline(0, 'k--');
    yline(0, 'k:');
    xlabel('Time (s)');
    ylabel('Beta Contra - Ipsi ERD (%)');
    title([cond_labels{c} ' hand']);
    grid on;
    hold off;
end

saveas(figure_beta_diff, fullfile(output_dir, 'ERD_Beta_ContraMinusIpsi.png'));
close(figure_beta_diff);

fprintf('[%s] STEP7 ERD analysis complete. Outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), output_dir);
