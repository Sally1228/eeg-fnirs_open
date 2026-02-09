% EEG-fNIRS block-level correlation for cate task (Category Verbal Fluency).
% Assumptions & Notes:
% - EEG epochs are time-locked to marker=1 task blocks.
% - EEG data file naming follows *_task-cate_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set
%   (fallback patterns are checked below).
% - fNIRS block averages are from STEP3_BlockAvarage and saved to *_task-cate_blockavg_table.csv.
% - Block order is derived from EEG event latency; fNIRS blockIndex is chronological within task.
% - If blocks are missing in either modality, use the intersection of block indices.
% - Correlations are Pearson; p-values require Statistics Toolbox (tcdf).
% - EEG ROI labels follow 10-5 naming; missing channels are skipped per ROI.
% - EEG ERD is computed for Delta/Theta/Alpha/Beta (task vs baseline).


% 脚本目的：计算 cate 任务下每个 block 的 EEG 频带 ERD（Delta/Theta/Alpha/Beta），并与同一 block 的 fNIRS HbO/HbR 均值做相关分析（按 ROI 和全体）。
% 关键假设
% - EEG epoch 对齐 marker=1（每个 epoch 是一个任务 block）。
% - fNIRS block 平均来自 STEP3 的 *_blockavg_table.csv。
% - EEG block 顺序按事件 latency 排序；fNIRS 用 blockIndex。
% - ROI 使用 cate 的四个 ROI：LeftFrontal / LeftTemporal / RightFrontal / RightTemporal。
% - ERD 计算为“任务功率 vs 基线功率”的百分比变化。
% 输入
% - EEG：eeg/cate/sub-xxx/ 下的 epoched .set 文件（脚本有多个文件名候选）。
% -fNIRS：sub-xxx_task-cate_blockavg_table.csv。
% 主要流程
% - 寻找共同被试：在 EEG 与 fNIRS 文件夹中取交集。
% - 逐被试处理：读取 EEG .set 与 fNIRS block 表。用 EEG 的时间轴生成 baseline 和 task 的时间窗口。根据 EEG 通道标签找到 ROI 通道索引。标记坏 trial（EEG.reject + boundary 事件）。用事件 latency 给每个 epoch 赋 blockIndex。
% - 逐 ROI 计算 EEG ERD：对每个 block 的 ROI 信号做 PSD。在 Delta/Theta/Alpha/Beta 频带分别计算 ERD%。
% - EEG 与 fNIRS block 对齐：fNIRS 按 blockIndex 提取 HbO/HbR blockMean。只保留 EEG 与 fNIRS 都存在的 block。
% - 相关分析：被试内：ROI‑level + All‑ROI，分别对每个频带与 HbO/HbR 做相关。组水平：所有被试合并后再做同样相关。
% 输出
% - cate_eeg_fnirs_block_features.csv包含每个 block 的 EEG ERD（4 个频带） + fNIRS HbO/HbR。
% - cate_eeg_fnirs_block_alignment.csv 对齐统计（EEG block 数、fNIRS block 数、匹配数）。
% - cate_eeg_fnirs_block_correlation.csv 相关系数 R、P 值等结果。
% - cate_eeg_fnirs_block_results.mat 保存表格和参数（包括 band_defs）。


%% Config
close all; clearvars; clc;

taskName = 'cate';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset patterns (epoched)
eegSetPatterns = {'%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'};

% fNIRS block average table
fnirsBlockTablePattern = '%s_task-%s_blockavg_table.csv';

% Single condition mapping
cond_label = 'task';
cond_marker = 1;

% EEG ROI labels
roi_names = {'LeftFrontal','LeftTemporal','RightFrontal','RightTemporal'};
roi_labels = {
    {'F1','F3','F5','AF3','FC1','FC3','FC5','C1','C3'}, ...
    {'FT7','T7','TP7','CP5','C5'}, ...
    {'F2','F4','F6','AF4','FC2','FC4','FC6','C2','C4'}, ...
    {'FT8','T8','TP8','CP6','C6'} ...
};

% EEG time windows (ms)
baseline_ms = [-2000 0];
task_window_ms = [0 60000];

% EEG bands (Delta/Theta/Alpha/Beta)
band_defs = {
    'Delta', [1 4];
    'Theta', [4 8];
    'Alpha', [8 13];
    'Beta',  [13 30]
};
nBand = size(band_defs, 1);
bandNames = band_defs(:,1);
bandVarNames = cell(nBand,1);
for bi = 1:nBand
    bandVarNames{bi} = sprintf('EEG_%sERD_Pct', bandNames{bi});
end
min_samples = 32;

% Correlation settings
minBlocksForCorr = 3;

%% Paths
% script_dir = '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg_fnirs_processing/cate_eeg_fnirs_processing';
script_dir = fileparts(mfilename('fullpath'));

repo_dir = fileparts(fileparts(script_dir));
eeg_root = fullfile(repo_dir, 'eeg', taskName);
fnirs_root = fullfile(repo_dir, 'fnirs', taskName);
out_root = fullfile(script_dir, 'output', 'STEP1_cate_block_correlation');
if exist(out_root, 'dir') ~= 7
    mkdir(out_root);
end

fprintf('[%s] Assumptions: marker=1 task blocks, block order by EEG latency, fNIRS blockIndex chronological.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] Assumptions: EEG ROI and fNIRS ROI = Left/Right Frontal/Temporal.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));

if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>

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
blockRows = cell(0, 5 + nBand + 2);
alignRows = cell(0, 6);

hasTcdf = exist('tcdf', 'file') == 2;
if ~hasTcdf
    fprintf('[%s] WARNING: tcdf not found; p-values will be NaN.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end

%% Main loop
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
    hasBlockLabel = ismember('blockLabel', fnirsT.Properties.VariableNames);

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
        fprintf('[%s] WARNING: No EEG ROI channels found for %s. Skipping subject.\n', ...
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

    %% Per-ROI processing
    cond_all_trials = find(epoch_bin == cond_marker);
    cond_good_trials = find(good_trials & epoch_bin == cond_marker);

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

    for r = 1:numel(roi_names)
        roi_idx = roi_idx_all{r};
        if isempty(roi_idx)
            continue;
        end
        roi_name = roi_names{r};

        % EEG block features for good trials
        eeg_block_idx = [];
        eeg_band_erd = [];

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

            band_erd = nan(1, nBand);
            for bi = 1:nBand
                band_range = band_defs{bi,2};
                base_idx = f_base >= band_range(1) & f_base <= band_range(2);
                task_idx = f_task >= band_range(1) & f_task <= band_range(2);
                base_pow = NaN;
                task_pow = NaN;
                if any(base_idx)
                    base_pow = trapz(f_base(base_idx), Pxx_base(base_idx));
                end
                if any(task_idx)
                    task_pow = trapz(f_task(task_idx), Pxx_task(task_idx));
                end
                if isfinite(base_pow) && base_pow > 0 && isfinite(task_pow)
                    band_erd(bi) = (task_pow / base_pow - 1) * 100;
                end
            end

            eeg_block_idx(end+1, 1) = block_idx; %#ok<AGROW>
            eeg_band_erd(end+1, 1:nBand) = band_erd; %#ok<AGROW>
        end

        % fNIRS block means
        condMask = strcmpi(fnirsT.condition, cond_label) & strcmpi(fnirsT.roi, roi_name);
        hboMask = condMask & strcmpi(fnirsT.chromophore, 'HbO');
        hbrMask = condMask & strcmpi(fnirsT.chromophore, 'HbR');

        block_idx_hbo = fnirsT.blockIndex(hboMask);
        block_idx_hbr = fnirsT.blockIndex(hbrMask);
        hbo_vals = fnirsT.blockMean(hboMask);
        hbr_vals = fnirsT.blockMean(hbrMask);

        hbo_labels = {};
        if hasBlockLabel
            hbo_labels = fnirsT.blockLabel(hboMask);
        end

        if ~isempty(block_idx_hbo)
            [hbo_unique, ~, hbo_grp] = unique(block_idx_hbo);
            hbo_vals = accumarray(hbo_grp, hbo_vals, [], @(x) mean(x, 'omitnan'));
            if hasBlockLabel
                hbo_label_unique = cell(numel(hbo_unique), 1);
                for ui = 1:numel(hbo_unique)
                    idxu = find(hbo_grp == ui, 1, 'first');
                    hbo_label_unique{ui} = hbo_labels{idxu};
                end
                hbo_labels = hbo_label_unique;
            end
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
        fnirs_labels = {};
        if hasBlockLabel && ~isempty(hbo_labels)
            fnirs_labels = hbo_labels(idx_hbo);
        end

        % Align EEG and fNIRS blocks
        [common_blocks, idx_fnirs, idx_eeg] = intersect(fnirs_block_idx, eeg_block_idx);

        alignRows(end+1, :) = {subName, roi_name, numel(cond_all_trials), numel(eeg_block_idx), ...
            numel(fnirs_block_idx), numel(common_blocks)}; %#ok<AGROW>

        for bi = 1:numel(common_blocks)
            block_label = '';
            if hasBlockLabel && ~isempty(fnirs_labels)
                block_label = fnirs_labels{idx_fnirs(bi)};
            end
            band_vals = eeg_band_erd(idx_eeg(bi), :);
            blockRows(end+1, :) = [{subName, cond_label, roi_name, common_blocks(bi), block_label}, ...
                num2cell(band_vals), {fnirs_hbo(idx_fnirs(bi)), fnirs_hbr(idx_fnirs(bi))}]; %#ok<AGROW>
        end
    end
end

%% Build tables
if isempty(blockRows)
    error('No matched blocks found. Check EEG/fNIRS inputs and ROI/marker settings.');
end

blockTable = cell2table(blockRows, 'VariableNames', ...
    [{'Subject','Condition','ROI','BlockIndex','BlockLabel'}, bandVarNames', ...
     {'FNIRS_HbO_BlockMean','FNIRS_HbR_BlockMean'}]);

alignmentTable = cell2table(alignRows, 'VariableNames', ...
    {'Subject','ROI','EEG_TotalBlocks','EEG_GoodBlocks','FNIRS_Blocks','MatchedBlocks'});

%% Correlation analysis
subjects = unique(blockTable.Subject, 'stable');
analysisNames = [roi_names {'All'}];
bandNames = band_defs(:,1);
eegCols = bandVarNames;
chromNames = {'HbO', 'HbR'};
fnirsCols = {'FNIRS_HbO_BlockMean', 'FNIRS_HbR_BlockMean'};

corrRows = cell(0, 8);
for si = 1:numel(subjects)
    subName = subjects{si};
    for ai = 1:numel(analysisNames)
        analysisName = analysisNames{ai};
        if strcmpi(analysisName, 'All')
            roiMask = true(height(blockTable), 1);
        else
            roiMask = strcmpi(blockTable.ROI, analysisName);
        end
        subMask = strcmp(blockTable.Subject, subName);
        mask = subMask & roiMask;
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
                    if strcmpi(analysisName, 'All')
                        roi_vals = blockTable.ROI(mask);
                        block_ids = blockTable.BlockIndex(mask);
                        roi_vals = roi_vals(valid);
                        block_ids = block_ids(valid);
                        detail_parts = cell(numel(block_ids), 1);
                        for di = 1:numel(block_ids)
                            detail_parts{di} = sprintf('%s:%d', roi_vals{di}, block_ids(di));
                        end
                        detail = strjoin(detail_parts, ';');
                    else
                        block_ids = blockTable.BlockIndex(mask);
                        block_ids = block_ids(valid);
                        detail = strjoin(arrayfun(@(v) sprintf('%d', v), block_ids(:)', 'UniformOutput', false), ';');
                    end
                end
                corrRows(end+1, :) = {subName, analysisName, bandNames{bi}, chromNames{ci}, n, r, p, detail}; %#ok<AGROW>
            end
        end
    end
end

% Group-level correlation across all subjects
for ai = 1:numel(analysisNames)
    analysisName = analysisNames{ai};
    if strcmpi(analysisName, 'All')
        roiMask = true(height(blockTable), 1);
    else
        roiMask = strcmpi(blockTable.ROI, analysisName);
    end
    idx_all = find(roiMask);
    for bi = 1:numel(bandNames)
        x = blockTable.(eegCols{bi})(roiMask);
        for ci = 1:numel(chromNames)
            y = blockTable.(fnirsCols{ci})(roiMask);
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
    {'Subject','ROI','EEGBand','FNIRSChrom','N','R','P','Detail'});

%% Save outputs
blockCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_features.csv', taskName));
alignCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_alignment.csv', taskName));
corrCsv = fullfile(out_root, sprintf('%s_eeg_fnirs_block_correlation.csv', taskName));

writetable(blockTable, blockCsv);
writetable(alignmentTable, alignCsv);
writetable(corrTable, corrCsv);

save(fullfile(out_root, sprintf('%s_eeg_fnirs_block_results.mat', taskName)), ...
    'blockTable', 'alignmentTable', 'corrTable', 'baseline_ms', 'task_window_ms', ...
    'band_defs', 'roi_names', 'roi_labels');

fprintf('[%s] Done. Outputs in %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), out_root);
