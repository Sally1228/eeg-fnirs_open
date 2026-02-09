% EEG-informed fNIRS GLM using single-trial P300 amplitude as a parametric modulator.
% Assumptions & Notes:
% - EEG datasets are epoched and binned with ERPLAB; common bin = 1, oddball bin = 2.
% - Artifact rejection marks are stored in EEG.reject and epochs are not removed.
% - fNIRS preproc MAT contains dc and stim; stim(1)=common, stim(2)=oddball.
% - EEG and fNIRS oddball/common events are aligned by order; if counts mismatch, time-based matching is attempted.
% - fNIRS trial validity uses tIncAuto (global); channel-wise tIncAutoCh is not used for trial filtering.
% - Oddball-common EEG regressors use the same weight scaling as all trials (common+oddball).


% 说明：用EEG单试次P300幅度（Pz）作为事件参数调制器，构建EEG‑informed fNIRS GLM，评估神经血管耦合（EEG振幅变化是否预测HbO/HbR变化）。
% 读取EEG与fNIRS的被试列表（交集），逐被试处理。
% EEG端：
% - 读取已分段/分bin的 .set 文件（优先带 _epoch_interp_ar 的版本）。
% - 按ERPLAB bin=1/2 找common/oddball试次，依据 EEG.reject 标记剔除坏试次（不直接删除）。
% - 在Pz计算单试次P300（300–500ms窗口，-200–0ms基线），得到每个试次的幅度。
% - 尝试从 EEG.urevent 取oddball/common事件时间（秒），用于与fNIRS对齐；若取不到则只能按顺序对齐。
% fNIRS端：
% - 读取 *_fnirs_preproc.mat，提取 dc、stim、tIncAuto、tIncAutoCh。
% - stim(1)=common, stim(2)=oddball，整理为 [onset, dur, amp] 并按时间排序。
% - 使用 tIncAuto 做trial级有效性筛选（默认事件窗口内≥80%时间为有效才保留）。
% EEG‑fNIRS对齐：
% - 若oddball/common数量一致，直接按顺序匹配。
% - 若数量不一致：用EEG事件时间与fNIRS oddball/common时间做时间匹配（容差0.3s），估计全局offset；匹配失败则回退到顺序匹配。
% 构建设计矩阵（3种分析）：
% - All auditory：所有刺激 → HRF卷积；所有刺激×EEG权重 → HRF卷积（核心）
% - Common-only：common → HRF卷积；common×EEG权重 → HRF卷积（核心）
% - Oddball-Common：oddball-common → HRF卷积；(oddball-common)×EEG权重 → HRF卷积（核心）
% - 加入截距 + 漂移项（默认DCT高通）
% GLM拟合：
% - 默认OLS，逐通道拟合；使用 tIncAuto/tIncAutoCh 做时间点掩码。
% - 输出每通道每回归量的β/t/p；可选ROI汇总。
% 关键配置（可调）
% - EEG相关：eegCommonBin=1, targetBin=2(oddball), eegChanLabel='Pz', p300WindowMs=[300 500], baselineWindowMs=[-200 0]
% - 事件对齐：alignToleranceSec=0.3, minTimeMatchEvents=5
% - fNIRS trial保留：fnirsTrialKeepFrac=0.8
% - HRF：gammaShape=6, gammaScale=1, hrfDurationSec=30
% - GLM漂移：driftModel='dct', dctCutoffSec=128
% - EEG权重缩放：eegWeightMode='zscore'
% - 可切换：useBlockAverageEEG=true 用block平均P300替代单试次
% 输出内容
% - 主结果：<sub-xxx>_task-audi_eegp300_glm_<analysis>.mat
% - 通道统计表：<sub-xxx>_task-audi_eegp300_glm_table_<analysis>.csv
% - 对齐表（检查匹配与权重）：<sub-xxx>_task-audi_eeg_fnirs_alignment.csv（oddball）
% - 对齐表（common）：<sub-xxx>_task-audi_eeg_fnirs_alignment_common.csv
% - ROI汇总（可选）：..._eegp300_glm_roi_table_<analysis>.csv
% - 图：HRF/回归量/权重分布、EEG‑调制β分布图（每个analysis）


%% Config
close all; clearvars; clc;

taskName = 'audi';
subjList = {}; % empty = auto (intersection of EEG and fNIRS subjects)

% EEG dataset pattern priority (first match wins)
eegSetPatterns = { ...
    '%s_task-%s_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set', ...
    };

% EEG oddball/common settings
eegCommonBin = 1;
targetBin = 2;
eegChanLabel = 'Pz';
eegChanFallback = 53;
p300WindowMs = [300 500];
baselineWindowMs = [-200 0];
p300Measure = 'peak'; % 'peak' | 'mean'
p300Polarity = 'positive'; % 'positive' | 'negative' | 'abs'
baselineMode = 'subtract'; % 'subtract' | 'none'

% Analysis settings
analysisNames = {'All','Common','OddballCommon'};
analysisTags = {'all','common','oddballcommon'};

% Alignment settings
alignToleranceSec = 0.3;
minTimeMatchEvents = 5;

% fNIRS event settings
stimIndexCommon = 1;
stimIndexOddball = 2;
stimDurationSec = 0.5;
fnirsTrialKeepFrac = 0.8; % require this fraction of tIncAuto==1 during event window
fnirsTrialPadSec = 0.0;   % extra padding before/after event window

% EEG weight settings
eegWeightMode = 'zscore'; % 'zscore' | 'demean' | 'raw'
useBlockAverageEEG = false;
blockGapSec = 10; % gap (sec) to detect block starts
blockTrials = 23; % total trials per block (common+oddball)
minValidTrialsPerBlock = 1;

% HRF settings (gamma)
hrfDurationSec = 30;
gammaShape = 6;
gammaScale = 1;

% GLM settings
glmMethod = 'ols'; % 'ols' | 'ar-irls'
driftModel = 'dct'; % 'none' | 'linear' | 'poly' | 'dct'
driftPolyOrder = 3;
dctCutoffSec = 128;

% ROI filtering (same as fnirs/audi/audi_fnirs_processing/STEP4_GLM.m)
useRoiFilter = true;
saveRoiSummary = true;
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

% Output settings
outputDirName = 'STEP2_p300_informed_fnirs';
saveFigures = true;
processOnlyFirst = false;

% Expected counts (optional QA)
expectedTotalTrials = 184;
expectedOddballTrials = 40;
expectedCommonTrials = 144;

%% Setup
% scriptDir = '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg_fnirs_processing/audi_eeg_fnirs_processing'
scriptDir = fileparts(mfilename('fullpath'));

rootDir = fileparts(fileparts(scriptDir));
eegRoot = fullfile(rootDir, 'eeg', taskName);
fnirsRoot = fullfile(rootDir, 'fnirs', taskName);
outRoot = fullfile(scriptDir, 'output', outputDirName);
if exist(outRoot, 'dir') ~= 7
    mkdir(outRoot);
end

fprintf('[%s] Assumptions: EEG bins common=%d oddball=%d; stim(1)=common stim(2)=oddball; tIncAuto trial QC; oddball-common uses all-trial EEG scaling.\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), eegCommonBin, targetBin);

if exist('eeglab', 'file') ~= 2
    error('EEGLAB not found on MATLAB path.');
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %#ok<ASGLU>

if strcmpi(glmMethod, 'ar-irls')
    if exist('ar_glm_final', 'file') ~= 2 || exist('robust_ar_fit', 'file') ~= 2
        error('AR-IRLS functions not found. Add Homer3 to path or use glmMethod="ols".');
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
    eegNames = {eegSubs.name};
    fnirsNames = {fnirsSubs.name};
    subjList = intersect(eegNames, fnirsNames, 'stable');
end
if isempty(subjList)
    error('No subjects found under %s and %s', eegRoot, fnirsRoot);
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
        fprintf('[%s] WARNING: EEG .set not found for %s. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    fnirsMat = fullfile(fnirsSubDir, sprintf('%s_task-%s_fnirs_preproc.mat', subName, taskName));
    if exist(fnirsMat, 'file') ~= 2
        fprintf('[%s] WARNING: fNIRS preproc MAT not found for %s. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    subOutDir = fullfile(outRoot, subName);
    if exist(subOutDir, 'dir') ~= 7
        mkdir(subOutDir);
    end

    %% Load EEG and compute single-trial P300
    EEG = pop_loadset('filename', eegSetFile, 'filepath', eegSubDir);

    if isempty(EEG.times) || EEG.trials < 1
        fprintf('[%s] WARNING: EEG data invalid for %s. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);
        continue;
    end

    chanIdx = [];
    if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
        labels = {EEG.chanlocs.labels};
        matchIdx = find(strcmpi(labels, eegChanLabel), 1);
        if ~isempty(matchIdx)
            chanIdx = matchIdx;
        end
    end
    if isempty(chanIdx)
        if eegChanFallback <= EEG.nbchan
            chanIdx = eegChanFallback;
            fprintf('[%s] WARNING: Channel %s not found; using index %d.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), eegChanLabel, chanIdx);
        else
            fprintf('[%s] WARNING: EEG channel not found and fallback invalid. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            continue;
        end
    end

    timeIdx = EEG.times >= p300WindowMs(1) & EEG.times <= p300WindowMs(2);
    baseIdx = EEG.times >= baselineWindowMs(1) & EEG.times <= baselineWindowMs(2);
    if ~any(timeIdx)
        fprintf('[%s] WARNING: P300 window outside EEG time range. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end

    % Build reject mask
    reject = false(1, EEG.trials);
    if isfield(EEG, 'reject') && ~isempty(EEG.reject)
        rejectFields = fieldnames(EEG.reject);
        for rf = 1:numel(rejectFields)
            val = EEG.reject.(rejectFields{rf});
            if isempty(val) || ~(islogical(val) || isnumeric(val))
                continue;
            end
            if isvector(val) && numel(val) == EEG.trials
                reject = reject | logical(val(:))';
            elseif ismatrix(val) && size(val, 2) == EEG.trials
                reject = reject | any(val, 1);
            elseif ismatrix(val) && size(val, 1) == EEG.trials
                reject = reject | any(val, 2)';
            end
        end
    end

    % Identify oddball/common epochs by bin
    oddEpochMask = false(1, EEG.trials);
    commonEpochMask = false(1, EEG.trials);
    if isfield(EEG, 'epoch') && ~isempty(EEG.epoch)
        for e = 1:EEG.trials
            if ~isfield(EEG.epoch(e), 'eventbini')
                continue;
            end
            eventBini = EEG.epoch(e).eventbini;
            biniVals = [];
            if iscell(eventBini)
                for bb = 1:numel(eventBini)
                    tmp = eventBini{bb};
                    if iscell(tmp)
                        tmp = [tmp{:}];
                    end
                    if ischar(tmp)
                        tmp = str2double(tmp);
                    end
                    if isempty(tmp)
                        continue;
                    end
                    biniVals = [biniVals, double(tmp(:)')]; %#ok<AGROW>
                end
            else
                biniVals = double(eventBini(:)');
            end
            if any(biniVals == targetBin)
                oddEpochMask(e) = true;
            end
            if any(biniVals == eegCommonBin)
                commonEpochMask(e) = true;
            end
        end
    end

    oddEpochIdx = find(oddEpochMask);
    commonEpochIdx = find(commonEpochMask);
    nOddEeg = numel(oddEpochIdx);
    nCommonEeg = numel(commonEpochIdx);
    if nOddEeg == 0 && nCommonEeg == 0
        fprintf('[%s] WARNING: No common/oddball epochs found. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    elseif nOddEeg == 0
        fprintf('[%s] WARNING: No oddball epochs found; oddball analyses will be skipped.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    elseif nCommonEeg == 0
        fprintf('[%s] WARNING: No common epochs found; common analyses will be skipped.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    end

    if expectedOddballTrials > 0 && nOddEeg ~= expectedOddballTrials
        fprintf('[%s] WARNING: EEG oddball trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nOddEeg, expectedOddballTrials);
    end
    if expectedCommonTrials > 0 && nCommonEeg ~= expectedCommonTrials
        fprintf('[%s] WARNING: EEG common trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nCommonEeg, expectedCommonTrials);
    end
    if expectedTotalTrials > 0 && (nOddEeg + nCommonEeg) ~= expectedTotalTrials
        fprintf('[%s] WARNING: EEG total trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), (nOddEeg + nCommonEeg), expectedTotalTrials);
    end

    eegOddballValid = ~reject(oddEpochIdx);
    eegCommonValid = ~reject(commonEpochIdx);

    dataChan = squeeze(EEG.data(chanIdx, :, :));
    if isvector(dataChan)
        dataChan = reshape(dataChan, [numel(EEG.times), 1]);
    end

    p300Vals = nan(nOddEeg, 1);
    if nOddEeg > 0
        for oi = 1:nOddEeg
            ep = oddEpochIdx(oi);
            if ~eegOddballValid(oi)
                continue;
            end
            trace = dataChan(:, ep);
            baseVal = 0;
            if strcmpi(baselineMode, 'subtract') && any(baseIdx)
                baseVal = mean(trace(baseIdx), 'omitnan');
            end

            winData = trace(timeIdx);
            if strcmpi(p300Measure, 'mean')
                peakVal = mean(winData, 'omitnan');
            else
                switch lower(p300Polarity)
                    case 'positive'
                        peakVal = max(winData);
                    case 'negative'
                        peakVal = min(winData);
                    otherwise
                        peakVal = max(abs(winData));
                end
            end
            p300Vals(oi) = peakVal - baseVal;
        end
    end

    p300ValsCommon = nan(nCommonEeg, 1);
    if nCommonEeg > 0
        for ci = 1:nCommonEeg
            ep = commonEpochIdx(ci);
            if ~eegCommonValid(ci)
                continue;
            end
            trace = dataChan(:, ep);
            baseVal = 0;
            if strcmpi(baselineMode, 'subtract') && any(baseIdx)
                baseVal = mean(trace(baseIdx), 'omitnan');
            end

            winData = trace(timeIdx);
            if strcmpi(p300Measure, 'mean')
                peakVal = mean(winData, 'omitnan');
            else
                switch lower(p300Polarity)
                    case 'positive'
                        peakVal = max(winData);
                    case 'negative'
                        peakVal = min(winData);
                    otherwise
                        peakVal = max(abs(winData));
                end
            end
            p300ValsCommon(ci) = peakVal - baseVal;
        end
    end

    % EEG oddball times (sec) using urevent latency if available
    eegOddballTimeSec = nan(nOddEeg, 1);
    if nOddEeg > 0 && isfield(EEG, 'urevent') && ~isempty(EEG.urevent)
        for oi = 1:nOddEeg
            ep = oddEpochIdx(oi);
            if ~isfield(EEG.epoch(ep), 'eventbini') || ~isfield(EEG.epoch(ep), 'eventurevent')
                continue;
            end
            eventBini = EEG.epoch(ep).eventbini;
            eventUrevent = EEG.epoch(ep).eventurevent;
            if ~iscell(eventBini)
                eventBini = num2cell(eventBini);
            end
            if ~iscell(eventUrevent)
                eventUrevent = num2cell(eventUrevent);
            end
            nEv = min(numel(eventBini), numel(eventUrevent));
            for ev = 1:nEv
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
                if any(double(tmp) == targetBin)
                    u = eventUrevent{ev};
                    if iscell(u)
                        u = [u{:}];
                    end
                    if isempty(u)
                        continue;
                    end
                    u = u(1);
                    if u >= 1 && u <= numel(EEG.urevent)
                        lat = EEG.urevent(u).latency;
                        if ~isempty(lat) && isfinite(lat)
                            eegOddballTimeSec(oi) = (double(lat) - 1) / EEG.srate;
                        end
                    end
                    break;
                end
            end
        end
    end

    % EEG common times (sec) using urevent latency if available
    eegCommonTimeSec = nan(nCommonEeg, 1);
    if nCommonEeg > 0 && isfield(EEG, 'urevent') && ~isempty(EEG.urevent)
        for ci = 1:nCommonEeg
            ep = commonEpochIdx(ci);
            if ~isfield(EEG.epoch(ep), 'eventbini') || ~isfield(EEG.epoch(ep), 'eventurevent')
                continue;
            end
            eventBini = EEG.epoch(ep).eventbini;
            eventUrevent = EEG.epoch(ep).eventurevent;
            if ~iscell(eventBini)
                eventBini = num2cell(eventBini);
            end
            if ~iscell(eventUrevent)
                eventUrevent = num2cell(eventUrevent);
            end
            nEv = min(numel(eventBini), numel(eventUrevent));
            for ev = 1:nEv
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
                if any(double(tmp) == eegCommonBin)
                    u = eventUrevent{ev};
                    if iscell(u)
                        u = [u{:}];
                    end
                    if isempty(u)
                        continue;
                    end
                    u = u(1);
                    if u >= 1 && u <= numel(EEG.urevent)
                        lat = EEG.urevent(u).latency;
                        if ~isempty(lat) && isfinite(lat)
                            eegCommonTimeSec(ci) = (double(lat) - 1) / EEG.srate;
                        end
                    end
                    break;
                end
            end
        end
    end

    %% Load fNIRS preproc
    S = load(fnirsMat, 'dc', 'stim', 'mlActAuto', 'tIncAuto', 'tIncAutoCh');
    if ~isfield(S, 'dc') || ~isfield(S, 'stim')
        fprintf('[%s] WARNING: fNIRS MAT missing dc/stim. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
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
        fprintf('[%s] WARNING: fNIRS time/data length mismatch. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end
    nT = numel(t);

    nStim = numel(stim);
    if nStim < max(stimIndexCommon, stimIndexOddball)
        fprintf('[%s] WARNING: stim has only %d entries. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), nStim);
        continue;
    end

    stimCommon = stim(stimIndexCommon).data;
    stimOddball = stim(stimIndexOddball).data;
    if isempty(stimCommon) || isempty(stimOddball)
        fprintf('[%s] WARNING: stim common/oddball empty. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end

    % Ensure Nx3 [onset dur amp]
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
    if expectedOddballTrials > 0 && nOddFnirs ~= expectedOddballTrials
        fprintf('[%s] WARNING: fNIRS oddball trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nOddFnirs, expectedOddballTrials);
    end
    if expectedCommonTrials > 0 && nCommonFnirs ~= expectedCommonTrials
        fprintf('[%s] WARNING: fNIRS common trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nCommonFnirs, expectedCommonTrials);
    end
    if expectedTotalTrials > 0 && (nCommonFnirs + nOddFnirs) ~= expectedTotalTrials
        fprintf('[%s] WARNING: fNIRS total trials=%d (expected %d).\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), (nCommonFnirs + nOddFnirs), expectedTotalTrials);
    end

    % fNIRS trial validity from tIncAuto
    fnirsOddballValid = true(nOddFnirs, 1);
    fnirsCommonValid = true(nCommonFnirs, 1);
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
            else
                tIncVec = [tIncVec; ones(nT - numel(tIncVec), 1)];
            end
        end
        for oi = 1:nOddFnirs
            onset = stimOddball(oi,1);
            dur = stimOddball(oi,2);
            if ~isfinite(dur) || dur <= 0
                dur = stimDurationSec;
            end
            winStart = onset - fnirsTrialPadSec;
            winStop = onset + dur + fnirsTrialPadSec;
            idx = t >= winStart & t <= winStop;
            if any(idx)
                keepFrac = mean(tIncVec(idx));
                fnirsOddballValid(oi) = keepFrac >= fnirsTrialKeepFrac;
            else
                fnirsOddballValid(oi) = false;
            end
        end
        for ci = 1:nCommonFnirs
            onset = stimCommon(ci,1);
            dur = stimCommon(ci,2);
            if ~isfinite(dur) || dur <= 0
                dur = stimDurationSec;
            end
            winStart = onset - fnirsTrialPadSec;
            winStop = onset + dur + fnirsTrialPadSec;
            idx = t >= winStart & t <= winStop;
            if any(idx)
                keepFrac = mean(tIncVec(idx));
                fnirsCommonValid(ci) = keepFrac >= fnirsTrialKeepFrac;
            else
                fnirsCommonValid(ci) = false;
            end
        end
    end

    %% Align EEG and fNIRS oddball events
    mapE = [];
    mapF = [];
    if nOddEeg == 0 || nOddFnirs == 0
        fprintf('[%s] WARNING: Oddball events missing; skipping oddball alignment.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    elseif nOddEeg == nOddFnirs
        mapE = (1:nOddEeg)';
        mapF = (1:nOddFnirs)';
        if any(isfinite(eegOddballTimeSec))
            dt = eegOddballTimeSec - stimOddball(:,1);
            dt = dt(isfinite(dt));
            if ~isempty(dt)
                dtMed = median(dt);
                dtMad = median(abs(dt - dtMed));
                fprintf('[%s] Alignment QA: median offset=%.3f s, MAD=%.3f s.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), dtMed, dtMad);
            end
        end
    else
        fprintf('[%s] WARNING: Oddball count mismatch (EEG=%d, fNIRS=%d). Using time-based matching.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nOddEeg, nOddFnirs);

        nFinite = sum(isfinite(eegOddballTimeSec));
        if nFinite >= minTimeMatchEvents
            eegTimesFinite = eegOddballTimeSec(isfinite(eegOddballTimeSec));
            nCandE = min(numel(eegTimesFinite), 5);
            nCandF = min(nOddFnirs, 5);
            offsets = [];
            for i = 1:nCandE
                for j = 1:nCandF
                    offsets(end+1) = eegTimesFinite(i) - stimOddball(j,1); %#ok<AGROW>
                end
            end
            bestOffset = offsets(1);
            bestMatch = -inf;
            for oc = 1:numel(offsets)
                offset = offsets(oc);
                eIdx = 1;
                fIdx = 1;
                count = 0;
                while eIdx <= nOddEeg && fIdx <= nOddFnirs
                    if ~isfinite(eegOddballTimeSec(eIdx))
                        eIdx = eIdx + 1;
                        continue;
                    end
                    dt = (eegOddballTimeSec(eIdx) - offset) - stimOddball(fIdx,1);
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
            fprintf('[%s] Time alignment offset=%.3f s (matches=%d)\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), bestOffset, bestMatch);

            eIdx = 1;
            fIdx = 1;
            while eIdx <= nOddEeg && fIdx <= nOddFnirs
                if ~isfinite(eegOddballTimeSec(eIdx))
                    eIdx = eIdx + 1;
                    continue;
                end
                dt = (eegOddballTimeSec(eIdx) - bestOffset) - stimOddball(fIdx,1);
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
            fprintf('[%s] WARNING: Not enough EEG event times; falling back to order.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            nMatch = min(nOddEeg, nOddFnirs);
            mapE = (1:nMatch)';
            mapF = (1:nMatch)';
        end
    end

    %% Align EEG and fNIRS common events
    mapECommon = [];
    mapFCommon = [];
    if nCommonEeg == 0 || nCommonFnirs == 0
        fprintf('[%s] WARNING: Common events missing; skipping common alignment.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    elseif nCommonEeg == nCommonFnirs
        mapECommon = (1:nCommonEeg)';
        mapFCommon = (1:nCommonFnirs)';
        if any(isfinite(eegCommonTimeSec))
            dt = eegCommonTimeSec - stimCommon(:,1);
            dt = dt(isfinite(dt));
            if ~isempty(dt)
                dtMed = median(dt);
                dtMad = median(abs(dt - dtMed));
                fprintf('[%s] Alignment QA (common): median offset=%.3f s, MAD=%.3f s.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), dtMed, dtMad);
            end
        end
    else
        fprintf('[%s] WARNING: Common count mismatch (EEG=%d, fNIRS=%d). Using time-based matching.\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), nCommonEeg, nCommonFnirs);

        nFinite = sum(isfinite(eegCommonTimeSec));
        if nFinite >= minTimeMatchEvents
            eegTimesFinite = eegCommonTimeSec(isfinite(eegCommonTimeSec));
            nCandE = min(numel(eegTimesFinite), 5);
            nCandF = min(nCommonFnirs, 5);
            offsets = [];
            for i = 1:nCandE
                for j = 1:nCandF
                    offsets(end+1) = eegTimesFinite(i) - stimCommon(j,1); %#ok<AGROW>
                end
            end
            bestOffset = offsets(1);
            bestMatch = -inf;
            for oc = 1:numel(offsets)
                offset = offsets(oc);
                eIdx = 1;
                fIdx = 1;
                count = 0;
                while eIdx <= nCommonEeg && fIdx <= nCommonFnirs
                    if ~isfinite(eegCommonTimeSec(eIdx))
                        eIdx = eIdx + 1;
                        continue;
                    end
                    dt = (eegCommonTimeSec(eIdx) - offset) - stimCommon(fIdx,1);
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
            fprintf('[%s] Time alignment (common) offset=%.3f s (matches=%d)\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), bestOffset, bestMatch);

            eIdx = 1;
            fIdx = 1;
            while eIdx <= nCommonEeg && fIdx <= nCommonFnirs
                if ~isfinite(eegCommonTimeSec(eIdx))
                    eIdx = eIdx + 1;
                    continue;
                end
                dt = (eegCommonTimeSec(eIdx) - bestOffset) - stimCommon(fIdx,1);
                if abs(dt) <= alignToleranceSec
                    mapECommon(end+1,1) = eIdx; %#ok<AGROW>
                    mapFCommon(end+1,1) = fIdx; %#ok<AGROW>
                    eIdx = eIdx + 1;
                    fIdx = fIdx + 1;
                elseif dt < -alignToleranceSec
                    eIdx = eIdx + 1;
                else
                    fIdx = fIdx + 1;
                end
            end
        else
            fprintf('[%s] WARNING: Not enough EEG event times (common); falling back to order.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            nMatch = min(nCommonEeg, nCommonFnirs);
            mapECommon = (1:nMatch)';
            mapFCommon = (1:nMatch)';
        end
    end

    %% Map EEG P300 to fNIRS oddball events
    eegWeightRaw = nan(nOddFnirs,1);
    eegEpochIndex = nan(nOddFnirs,1);
    eegOddballIndex = nan(nOddFnirs,1);
    eegRejected = nan(nOddFnirs,1);
    eegTimeForFnirs = nan(nOddFnirs,1);

    for mi = 1:numel(mapE)
        eIdx = mapE(mi);
        fIdx = mapF(mi);
        eegOddballIndex(fIdx) = eIdx;
        eegEpochIndex(fIdx) = oddEpochIdx(eIdx);
        eegRejected(fIdx) = double(~eegOddballValid(eIdx));
        eegTimeForFnirs(fIdx) = eegOddballTimeSec(eIdx);
        if eegOddballValid(eIdx) && fnirsOddballValid(fIdx)
            eegWeightRaw(fIdx) = p300Vals(eIdx);
        end
    end

    %% Map EEG P300 to fNIRS common events
    eegWeightRawCommon = nan(nCommonFnirs,1);
    eegEpochIndexCommon = nan(nCommonFnirs,1);
    eegCommonIndex = nan(nCommonFnirs,1);
    eegRejectedCommon = nan(nCommonFnirs,1);
    eegTimeForFnirsCommon = nan(nCommonFnirs,1);

    for mi = 1:numel(mapECommon)
        eIdx = mapECommon(mi);
        fIdx = mapFCommon(mi);
        eegCommonIndex(fIdx) = eIdx;
        eegEpochIndexCommon(fIdx) = commonEpochIdx(eIdx);
        eegRejectedCommon(fIdx) = double(~eegCommonValid(eIdx));
        eegTimeForFnirsCommon(fIdx) = eegCommonTimeSec(eIdx);
        if eegCommonValid(eIdx) && fnirsCommonValid(fIdx)
            eegWeightRawCommon(fIdx) = p300ValsCommon(eIdx);
        end
    end

    if useBlockAverageEEG
        onsetsAll = sort([stimCommon(:,1); stimOddball(:,1)]);
        gapIdx = find(diff(onsetsAll) > blockGapSec);
        blockStartTimes = onsetsAll([1; gapIdx + 1]);
        blockIndexOddball = zeros(nOddFnirs,1);
        blockIndexCommon = zeros(nCommonFnirs,1);
        if numel(blockStartTimes) >= 2
            for oi = 1:nOddFnirs
                blockIndexOddball(oi) = sum(stimOddball(oi,1) >= blockStartTimes);
            end
            for ci = 1:nCommonFnirs
                blockIndexCommon(ci) = sum(stimCommon(ci,1) >= blockStartTimes);
            end
        else
            mergedOnsets = [stimCommon(:,1); stimOddball(:,1)];
            mergedIsOddball = [zeros(nCommonFnirs,1); ones(nOddFnirs,1)];
            [~, sortIdx] = sort(mergedOnsets);
            mergedIsOddball = mergedIsOddball(sortIdx);
            trialIdxOddball = find(mergedIsOddball == 1);
            trialIdxCommon = find(mergedIsOddball == 0);
            blockIndexOddball = ceil(trialIdxOddball / blockTrials);
            blockIndexCommon = ceil(trialIdxCommon / blockTrials);
        end
        nBlocks = max([blockIndexOddball(:); blockIndexCommon(:); 0]);
        blockMeanOdd = nan(nBlocks,1);
        blockMeanCommon = nan(nBlocks,1);
        for bi = 1:nBlocks
            maskOdd = blockIndexOddball == bi & isfinite(eegWeightRaw);
            if sum(maskOdd) >= minValidTrialsPerBlock
                blockMeanOdd(bi) = mean(eegWeightRaw(maskOdd), 'omitnan');
            end
            maskCommon = blockIndexCommon == bi & isfinite(eegWeightRawCommon);
            if sum(maskCommon) >= minValidTrialsPerBlock
                blockMeanCommon(bi) = mean(eegWeightRawCommon(maskCommon), 'omitnan');
            end
        end
        for oi = 1:nOddFnirs
            bi = blockIndexOddball(oi);
            if bi >= 1 && bi <= numel(blockMeanOdd) && isfinite(blockMeanOdd(bi))
                eegWeightRaw(oi) = blockMeanOdd(bi);
            else
                eegWeightRaw(oi) = NaN;
            end
        end
        for ci = 1:nCommonFnirs
            bi = blockIndexCommon(ci);
            if bi >= 1 && bi <= numel(blockMeanCommon) && isfinite(blockMeanCommon(bi))
                eegWeightRawCommon(ci) = blockMeanCommon(bi);
            else
                eegWeightRawCommon(ci) = NaN;
            end
        end
    end

    validWeightOdd = isfinite(eegWeightRaw);
    validWeightCommon = isfinite(eegWeightRawCommon);

    eegWeightUsed = eegWeightRaw;
    weightMean = NaN;
    weightStd = NaN;
    if any(validWeightOdd)
        weightMean = mean(eegWeightUsed(validWeightOdd));
        weightStd = std(eegWeightUsed(validWeightOdd));
        switch lower(eegWeightMode)
            case 'zscore'
                if weightStd <= 0 || ~isfinite(weightStd)
                    weightStd = 1;
                end
                eegWeightUsed(validWeightOdd) = (eegWeightUsed(validWeightOdd) - weightMean) / weightStd;
            case 'demean'
                eegWeightUsed(validWeightOdd) = eegWeightUsed(validWeightOdd) - weightMean;
            case 'raw'
                % no scaling
            otherwise
                fprintf('[%s] WARNING: Unknown eegWeightMode; using raw.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
    end
    eegWeightUsed(~validWeightOdd) = 0;

    eegWeightUsedCommon = eegWeightRawCommon;
    weightMeanCommon = NaN;
    weightStdCommon = NaN;
    if any(validWeightCommon)
        weightMeanCommon = mean(eegWeightUsedCommon(validWeightCommon));
        weightStdCommon = std(eegWeightUsedCommon(validWeightCommon));
        switch lower(eegWeightMode)
            case 'zscore'
                if weightStdCommon <= 0 || ~isfinite(weightStdCommon)
                    weightStdCommon = 1;
                end
                eegWeightUsedCommon(validWeightCommon) = (eegWeightUsedCommon(validWeightCommon) - weightMeanCommon) / weightStdCommon;
            case 'demean'
                eegWeightUsedCommon(validWeightCommon) = eegWeightUsedCommon(validWeightCommon) - weightMeanCommon;
            case 'raw'
                % no scaling
            otherwise
                fprintf('[%s] WARNING: Unknown eegWeightMode; using raw.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
    end
    eegWeightUsedCommon(~validWeightCommon) = 0;

    eegWeightUsedAllOdd = eegWeightRaw;
    eegWeightUsedAllCommon = eegWeightRawCommon;
    weightMeanAll = NaN;
    weightStdAll = NaN;
    allWeights = [eegWeightRawCommon(:); eegWeightRaw(:)];
    validAll = isfinite(allWeights);
    if any(validAll)
        weightMeanAll = mean(allWeights(validAll));
        weightStdAll = std(allWeights(validAll));
        switch lower(eegWeightMode)
            case 'zscore'
                if weightStdAll <= 0 || ~isfinite(weightStdAll)
                    weightStdAll = 1;
                end
                if any(validWeightCommon)
                    eegWeightUsedAllCommon(validWeightCommon) = (eegWeightUsedAllCommon(validWeightCommon) - weightMeanAll) / weightStdAll;
                end
                if any(validWeightOdd)
                    eegWeightUsedAllOdd(validWeightOdd) = (eegWeightUsedAllOdd(validWeightOdd) - weightMeanAll) / weightStdAll;
                end
            case 'demean'
                if any(validWeightCommon)
                    eegWeightUsedAllCommon(validWeightCommon) = eegWeightUsedAllCommon(validWeightCommon) - weightMeanAll;
                end
                if any(validWeightOdd)
                    eegWeightUsedAllOdd(validWeightOdd) = eegWeightUsedAllOdd(validWeightOdd) - weightMeanAll;
                end
            case 'raw'
                % no scaling
            otherwise
                fprintf('[%s] WARNING: Unknown eegWeightMode; using raw.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
    end
    eegWeightUsedAllCommon(~validWeightCommon) = 0;
    eegWeightUsedAllOdd(~validWeightOdd) = 0;

    fprintf('[%s] EEG oddball valid=%d/%d; EEG common valid=%d/%d; fNIRS oddball valid=%d/%d; fNIRS common valid=%d/%d; matched odd=%d; matched common=%d\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), sum(eegOddballValid), nOddEeg, sum(eegCommonValid), nCommonEeg, ...
        sum(fnirsOddballValid), nOddFnirs, sum(fnirsCommonValid), nCommonFnirs, numel(mapE), numel(mapECommon));

    %% Build HRF regressors on fNIRS time base
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        fprintf('[%s] WARNING: Invalid fNIRS time vector. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end
    tHrf = (0:dt:hrfDurationSec)';
    hrf = (tHrf.^(gammaShape-1)) .* exp(-tHrf/gammaScale) ./ (gamma(gammaShape) * (gammaScale^gammaShape));
    if max(hrf) > 0
        hrf = hrf ./ max(hrf);
    end

    uCommon = zeros(nT,1);
    for ei = 1:nCommonFnirs
        onset = stimCommon(ei,1);
        dur = stimCommon(ei,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        uCommon(idx) = uCommon(idx) + 1;
    end
    xCommon = conv(uCommon, hrf);
    xCommon = xCommon(1:nT);

    uOddball = zeros(nT,1);
    for ei = 1:nOddFnirs
        onset = stimOddball(ei,1);
        dur = stimOddball(ei,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        uOddball(idx) = uOddball(idx) + 1;
    end
    xOddball = conv(uOddball, hrf);
    xOddball = xOddball(1:nT);

    xAll = xCommon + xOddball;
    xDiff = xOddball - xCommon;

    uCommonEegAll = zeros(nT,1);
    for ei = 1:nCommonFnirs
        onset = stimCommon(ei,1);
        dur = stimCommon(ei,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        if eegWeightUsedAllCommon(ei) ~= 0
            uCommonEegAll(idx) = uCommonEegAll(idx) + eegWeightUsedAllCommon(ei);
        end
    end
    xCommonEegAll = conv(uCommonEegAll, hrf);
    xCommonEegAll = xCommonEegAll(1:nT);

    uOddballEegAll = zeros(nT,1);
    for ei = 1:nOddFnirs
        onset = stimOddball(ei,1);
        dur = stimOddball(ei,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        if eegWeightUsedAllOdd(ei) ~= 0
            uOddballEegAll(idx) = uOddballEegAll(idx) + eegWeightUsedAllOdd(ei);
        end
    end
    xOddballEegAll = conv(uOddballEegAll, hrf);
    xOddballEegAll = xOddballEegAll(1:nT);

    uCommonEegCommon = zeros(nT,1);
    for ei = 1:nCommonFnirs
        onset = stimCommon(ei,1);
        dur = stimCommon(ei,2);
        if ~isfinite(dur) || dur <= 0
            dur = stimDurationSec;
        end
        idx = t >= onset & t < (onset + dur);
        if eegWeightUsedCommon(ei) ~= 0
            uCommonEegCommon(idx) = uCommonEegCommon(idx) + eegWeightUsedCommon(ei);
        end
    end
    xCommonEegCommon = conv(uCommonEegCommon, hrf);
    xCommonEegCommon = xCommonEegCommon(1:nT);

    xAllEeg = xCommonEegAll + xOddballEegAll;
    xDiffEeg = xOddballEegAll - xCommonEegAll;

    %% ROI selection and QC masks
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
            fprintf('[%s] WARNING: No channels matched ROI filter. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            continue;
        end
    else
        roiLabel(:) = "ALL";
    end

    qcMaskAll = true(nMeasAll,1);
    if isfield(S, 'mlActAuto') && ~isempty(S.mlActAuto)
        if iscell(S.mlActAuto)
            mlActVec = S.mlActAuto{1};
        else
            mlActVec = S.mlActAuto;
        end
        mlActVec = mlActVec(:) ~= 0;
        if numel(mlActVec) == nMeasAll
            qcMaskAll = mlActVec;
        end
    end

    glmMask = roiMask & qcMaskAll;
    if ~any(glmMask)
        fprintf('[%s] WARNING: No channels after ROI/QC mask. Skipping.\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        continue;
    end

    y = y(:, glmMask);
    src = src(glmMask);
    det = det(glmMask);
    chrom = chrom(glmMask);
    roiLabel = roiLabel(glmMask);
    nMeas = numel(src);

    % Time masks (tIncAuto and tIncAutoCh)
    tIncCh = [];
    if isfield(S, 'tIncAutoCh') && ~isempty(S.tIncAutoCh)
        if iscell(S.tIncAutoCh)
            tIncCh = S.tIncAutoCh{1};
        else
            tIncCh = S.tIncAutoCh;
        end
        if size(tIncCh,1) ~= nT && size(tIncCh,2) == nT
            tIncCh = tIncCh';
        end
        if size(tIncCh,1) == nT
            if size(tIncCh,2) == nMeasAll
                tIncCh = double(tIncCh > 0);
                tIncCh = tIncCh(:, glmMask);
            else
                tIncCh = [];
            end
        else
            tIncCh = [];
        end
    end

    %% Alignment tables
    if nOddFnirs > 0
        alignmentTable = table();
        alignmentTable.oddIndexFnirs = (1:nOddFnirs)';
        alignmentTable.oddIndexEeg = eegOddballIndex;
        alignmentTable.eegEpochIndex = eegEpochIndex;
        alignmentTable.eegRejected = eegRejected;
        alignmentTable.fnirsBadTime = double(~fnirsOddballValid);
        alignmentTable.eegTimeSec = eegTimeForFnirs;
        alignmentTable.fnirsTimeSec = stimOddball(:,1);
        alignmentTable.eegP300 = eegWeightRaw;
        alignmentTable.eegWeightUsed = eegWeightUsed;
        alignCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_fnirs_alignment.csv', subName, taskName));
        writetable(alignmentTable, alignCsv);
    end

    if nCommonFnirs > 0
        alignmentTableCommon = table();
        alignmentTableCommon.commonIndexFnirs = (1:nCommonFnirs)';
        alignmentTableCommon.commonIndexEeg = eegCommonIndex;
        alignmentTableCommon.eegEpochIndex = eegEpochIndexCommon;
        alignmentTableCommon.eegRejected = eegRejectedCommon;
        alignmentTableCommon.fnirsBadTime = double(~fnirsCommonValid);
        alignmentTableCommon.eegTimeSec = eegTimeForFnirsCommon;
        alignmentTableCommon.fnirsTimeSec = stimCommon(:,1);
        alignmentTableCommon.eegP300 = eegWeightRawCommon;
        alignmentTableCommon.eegWeightUsed = eegWeightUsedCommon;
        alignCommonCsv = fullfile(subOutDir, sprintf('%s_task-%s_eeg_fnirs_alignment_common.csv', subName, taskName));
        writetable(alignmentTableCommon, alignCommonCsv);
    end

    %% Analyses
    if numel(analysisNames) ~= numel(analysisTags)
        error('analysisNames and analysisTags must have the same length.');
    end

    for ai = 1:numel(analysisNames)
        analysisName = analysisNames{ai};
        analysisTag = analysisTags{ai};
        xBase = [];
        xEeg = [];
        baseName = '';
        eegName = '';
        analysisWeightUsed = [];
        analysisWeightMean = NaN;
        analysisWeightStd = NaN;
        analysisTitle = analysisName;
        baseLabel = '';
        eegLabel = '';

        if strcmpi(analysisName, 'All') && (nCommonFnirs + nOddFnirs) == 0
            fprintf('[%s] WARNING: No auditory events; skipping %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
            continue;
        elseif strcmpi(analysisName, 'Common') && nCommonFnirs == 0
            fprintf('[%s] WARNING: No common events; skipping %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
            continue;
        elseif strcmpi(analysisName, 'OddballCommon') && (nCommonFnirs == 0 || nOddFnirs == 0)
            fprintf('[%s] WARNING: Missing oddball/common events; skipping %s.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
            continue;
        end

        switch lower(analysisName)
            case 'all'
                xBase = xAll;
                xEeg = xAllEeg;
                baseName = 'all';
                eegName = 'all_eeg';
                analysisWeightUsed = [eegWeightUsedAllCommon(validWeightCommon); eegWeightUsedAllOdd(validWeightOdd)];
                analysisWeightMean = weightMeanAll;
                analysisWeightStd = weightStdAll;
                baseLabel = 'All';
                eegLabel = 'All-EEG';
            case 'common'
                xBase = xCommon;
                xEeg = xCommonEegCommon;
                baseName = 'common';
                eegName = 'common_eeg';
                analysisWeightUsed = eegWeightUsedCommon(validWeightCommon);
                analysisWeightMean = weightMeanCommon;
                analysisWeightStd = weightStdCommon;
                baseLabel = 'Common';
                eegLabel = 'Common-EEG';
            case 'oddballcommon'
                xBase = xDiff;
                xEeg = xDiffEeg;
                baseName = 'oddball_common';
                eegName = 'oddball_common_eeg';
                analysisWeightUsed = [eegWeightUsedAllCommon(validWeightCommon); eegWeightUsedAllOdd(validWeightOdd)];
                analysisWeightMean = weightMeanAll;
                analysisWeightStd = weightStdAll;
                analysisTitle = 'Oddball-Common';
                baseLabel = 'Oddball-Common';
                eegLabel = 'Oddball-Common-EEG';
            otherwise
                fprintf('[%s] WARNING: Unknown analysisName=%s; skipping.\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
                continue;
        end

        if norm(xBase) <= 0
            fprintf('[%s] WARNING: %s base regressor is all zeros; skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
            continue;
        end
        if norm(xEeg) <= 0
            fprintf('[%s] WARNING: %s EEG regressor is all zeros; skipping.\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), analysisName);
            continue;
        end

        analysisWeightUsed = analysisWeightUsed(isfinite(analysisWeightUsed));

        % Design matrix
        X = [xBase xEeg];
        regNames = {baseName, eegName};

        % Intercept
        X = [X ones(nT,1)];
        regNames = [regNames {'intercept'}];

        % Drift terms
        switch lower(strtrim(driftModel))
            case 'none'
                % no drift
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
            tMask = tMask & all(isfinite(X), 2) & isfinite(y(:,mi));

            nKeep = sum(tMask);
            if nKeep <= p
                continue;
            end

            Xi = X(tMask, :);
            yi = y(tMask, mi);

            switch lower(glmMethod)
                case 'ar-irls'
                    Pmax = max(1, round(0.5 / dt));
                    [~, beta_i, tstat_i, pval_i, ~, CovB_i, dfe_i] = ar_glm_final(yi, Xi, Pmax);
                    CovB_i = squeeze(CovB_i);
                otherwise
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
                rowCell(end+1,:) = {subName, taskName, mi, src(mi), det(mi), char(chrom(mi)), char(roiLabel(mi)), ...
                    regNames{ri}, beta(ri,mi), seReg(ri,mi), tReg(ri,mi), pReg(ri,mi)}; %#ok<AGROW>
            end
        end

        T = cell2table(rowCell, 'VariableNames', ...
            {'subjId','task','measIndex','src','det','chromophore','roi','term','beta','SE','t','p'});

        outMat = fullfile(subOutDir, sprintf('%s_task-%s_eegp300_glm_%s.mat', subName, taskName, analysisTag));
        outCsv = fullfile(subOutDir, sprintf('%s_task-%s_eegp300_glm_table_%s.csv', subName, taskName, analysisTag));
        save(outMat, 'analysisName', 'analysisTag', 'beta', 'seReg', 'tReg', 'pReg', 'dof', 'regNames', 'X', ...
            'analysisWeightUsed', 'analysisWeightMean', 'analysisWeightStd', ...
            'eegWeightRaw', 'eegWeightUsed', 'weightMean', 'weightStd', ...
            'eegWeightRawCommon', 'eegWeightUsedCommon', 'weightMeanCommon', 'weightStdCommon', ...
            'eegWeightUsedAllOdd', 'eegWeightUsedAllCommon', 'weightMeanAll', 'weightStdAll', ...
            'chrom', 'src', 'det', 'roiLabel', 'roiMask', 'glmMask', 'qcMaskAll', ...
            'p300Vals', 'p300ValsCommon', 'eegOddballValid', 'eegCommonValid', ...
            'fnirsOddballValid', 'fnirsCommonValid', '-v7.3');
        writetable(T, outCsv);

        %% ROI summary (optional)
        if saveRoiSummary
            roiList = {};
            if useRoiFilter
                roiList = {roiDefs.name};
            else
                roiList = {'ALL'};
            end
            chromOrder = ["HbO","HbR","HbT"];
            chromList = chromOrder(ismember(chromOrder, unique(chrom)));
            termList = regNames;
            roiRows = {};
            for ri = 1:numel(roiList)
                roiName = roiList{ri};
                for ci = 1:numel(chromList)
                    chromName = chromList(ci);
                    chanIdx = (roiLabel == roiName) & (chrom == chromName);
                    for ti = 1:numel(termList)
                        termName = termList{ti};
                        regIdx = find(strcmp(regNames, termName), 1);
                        if isempty(regIdx)
                            continue;
                        end
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
            outRoiCsv = fullfile(subOutDir, sprintf('%s_task-%s_eegp300_glm_roi_table_%s.csv', subName, taskName, analysisTag));
            writetable(roiSummaryTable, outRoiCsv);
            save(outMat, 'roiSummaryTable', '-append');
        end

        %% Figures
        if saveFigures
            fig1 = figure('Visible','off');
            subplot(3,1,1);
            plot(tHrf, hrf, 'k', 'LineWidth', 1.5);
            xlabel('Time (s)'); ylabel('HRF');
            title('HRF (gamma)', 'Interpreter','none');

            subplot(3,1,2);
            plot(t, xBase, 'b'); hold on;
            plot(t, xEeg, 'g');
            legend({baseLabel, eegLabel}, 'Location','best');
            xlabel('Time (s)'); ylabel('Regressor');
            title(sprintf('Convolved regressors (%s)', analysisTitle), 'Interpreter','none');

            subplot(3,1,3);
            if ~isempty(analysisWeightUsed)
                histogram(analysisWeightUsed);
            else
                histogram(0);
            end
            xlabel('EEG weight (scaled)'); ylabel('Count');
            title('EEG weight distribution', 'Interpreter','none');

            figPath1 = fullfile(subOutDir, sprintf('%s_task-%s_eegp300_regressors_%s.png', subName, taskName, analysisTag));
            saveas(fig1, figPath1);
            close(fig1);

            eegRegIdx = find(strcmp(regNames, eegName), 1);
            if ~isempty(eegRegIdx)
                fig2 = figure('Visible','off');
                hboIdx = chrom == "HbO";
                hbrIdx = chrom == "HbR";
                hbtIdx = chrom == "HbT";
                dataPlot = [];
                grp = {};
                if any(hboIdx)
                    dataPlot = [dataPlot; beta(eegRegIdx, hboIdx)'];
                    grp = [grp; repmat({'HbO'}, sum(hboIdx), 1)];
                end
                if any(hbrIdx)
                    dataPlot = [dataPlot; beta(eegRegIdx, hbrIdx)'];
                    grp = [grp; repmat({'HbR'}, sum(hbrIdx), 1)];
                end
                if any(hbtIdx)
                    dataPlot = [dataPlot; beta(eegRegIdx, hbtIdx)'];
                    grp = [grp; repmat({'HbT'}, sum(hbtIdx), 1)];
                end
                if ~isempty(dataPlot)
                    boxplot(dataPlot, grp);
                else
                    plot(beta(eegRegIdx, :), 'k.');
                end
                ylabel(sprintf('Beta (%s)', eegName));
                title(sprintf('%s %s %s EEG-modulated beta', subName, taskName, analysisTitle), 'Interpreter','none');
                figPath2 = fullfile(subOutDir, sprintf('%s_task-%s_eegp300_beta_%s.png', subName, taskName, analysisTag));
                saveas(fig2, figPath2);
                close(fig2);
            end
        end
    end

    fprintf('[%s] Done: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), subName);

    if processOnlyFirst
        break;
    end
end
