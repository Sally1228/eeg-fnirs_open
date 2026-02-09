%==========================================================================
% 本脚本流程（概述）：
% 1) 输入/路径：读取各被试 ERP 文件列表与输出目录。
% 2) 预处理步骤：无额外预处理，直接基于已平均的 ERP 进行测量。
% 3) 核心计算：按指定时间窗与 ROI 通道测量 FN400 与 LPC 的幅度/潜伏期指标。
% 4) 输出：生成各项测量结果文本文件，保存到 outputs 目录下。
%==========================================================================
%Script #13
%Operates on individual subject data
%Uses the output from Script #8: Average_ERPs.m
%This script uses the individual subject averaged ERP waveforms from Script #7, and measures the mean amplitude, peak amplitude, peak latency, 50% area latency, and onset latency (50% peak latency) 
%during the time window of the component, and saves a separate text file for each measurement in the ERP Measurements folder. 
%Note that based on their respective susceptibility to high frequency noise, some measurements (e.g., mean amplitude, 50% area latency) are calculated on the averaged ERP waveforms without a low-pass filter applied, 
%whereas other measurements (e.g., peak amplitude, peak latency, and onset latency) are calculated from waveforms that have been low-pass filtered. 

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject ERP files are stored in the sibling data directory (eeg/recog/sub-XXX).
% - ERP channel labels include ROI names used for FN400/LPC.

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'outputs');
data_dir = fileparts(project_dir);

fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] Assumption: ERP channel labels include ROI names for FN400/LPC\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('[%s] STEP13 Measure ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to measure component amplitudes and latencies from (i.e., subjects that were not excluded due to excessive artifacts), based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-008'}; 
task ='recog';

%*************************************************************************************************************************************
%% FN400 
%Set measurement time window for measuring mean amplitude, peak amplitude, peak latency, and 50% area latency in milliseconds (e.g., 300 to 600 ms)
timewindow_fn400 = [300 500];
%Set measurement time window for measuring onset latency (50% peak latency) in milliseconds (e.g., 200 to 600 ms)
timewindow_onsetlat_fn400 = [100 500];
%Set EEG channel(s) to measure the components (ROI labels)
chan_labels_fn400 = {'Fz','FCz','Cz','F1','F2','FC1','FC2','C1','C2'};

%% LPC
%Set measurement time window for measuring mean amplitude, peak amplitude, peak latency, and 50% area latency in milliseconds
timewindow_lpc = [500 800];
%Set measurement time window for measuring onset latency (50% peak latency) in milliseconds
timewindow_onsetlat_lpc = [300 800];
%Set EEG channel(s) to measure the components (ROI labels)
chan_labels_lpc = {'P1','Pz','P2','CP1','CP2'};
%ROI average output
do_roi_average = true;
roi_output_dir_name = 'ROI_ERPsets';
roi_label_fn400 = 'ROI_FN400';
roi_label_lpc = 'ROI_LPC';
%Set difference wave bin(s) for measurement
diffbin = [3];  
%Set parent wave bins for measurement
parentbins = [1 2];  
%Set baseline correction period for measurement
baselinecorr = [-200 0]; 
%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%*************************************************************************************************************************************

%% Resolve channel indices from labels
labels_all = {};
labels_all_upper = {};
chanlocs_found = false;
for i = 1:length(SUB)
    Subject_Path = fullfile(data_dir, SUB{i});
    erppath = fullfile(Subject_Path, [SUB{i} '_recog_erp_ar_diff_waves.erp']);
    if exist(erppath, 'file')
        ERP = pop_loaderp('filename', [SUB{i} '_recog_erp_ar_diff_waves.erp'], 'filepath', Subject_Path);
        if isfield(ERP, 'chanlocs') && ~isempty(ERP.chanlocs)
            labels_all = {ERP.chanlocs.labels};
            labels_all = cellfun(@strtrim, labels_all, 'UniformOutput', false);
            labels_all_upper = cellfun(@upper, labels_all, 'UniformOutput', false);
            chanlocs_found = true;
            break;
        end
    end
end
if ~chanlocs_found
    error('No ERP file found to map channel labels.');
end

chan_labels_fn400_upper = cellfun(@upper, chan_labels_fn400, 'UniformOutput', false);
[~, idx_fn400] = ismember(chan_labels_fn400_upper, labels_all_upper);
missing_fn400 = chan_labels_fn400(idx_fn400 == 0);
if ~isempty(missing_fn400)
    warning('FN400 missing channels: %s', strjoin(missing_fn400, ', '));
end
chan_fn400 = idx_fn400(idx_fn400 > 0);
if isempty(chan_fn400)
    error('FN400 channel list is empty after lookup.');
end
fprintf('[%s] FN400 channels: %s -> %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    strjoin(chan_labels_fn400, ','), mat2str(chan_fn400));

chan_labels_lpc_upper = cellfun(@upper, chan_labels_lpc, 'UniformOutput', false);
[~, idx_lpc] = ismember(chan_labels_lpc_upper, labels_all_upper);
missing_lpc = chan_labels_lpc(idx_lpc == 0);
if ~isempty(missing_lpc)
    warning('LPC missing channels: %s', strjoin(missing_lpc, ', '));
end
chan_lpc = idx_lpc(idx_lpc > 0);
if isempty(chan_lpc)
    error('LPC channel list is empty after lookup.');
end
fprintf('[%s] LPC channels: %s -> %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    strjoin(chan_labels_lpc, ','), mat2str(chan_lpc));

%*************************************************************************************************************************************

%% Build ROI-averaged ERPsets (append ROI channels)
if do_roi_average
    roi_output_dir = fullfile(output_dir, roi_output_dir_name);
    if ~exist(roi_output_dir, 'dir')
        mkdir(roi_output_dir);
    end

    roi_chan_fn400 = numel(labels_all) + 1;
    roi_chan_lpc = numel(labels_all) + 2;

    roi_input_suffixes = { ...
        '_recog_erp_ar_diff_waves.erp', ...
        '_recog_erp_ar_diff_waves_lpfilt.erp' ...
        };
    roi_output_suffixes = { ...
        '_recog_erp_ar_diff_waves_roi.erp', ...
        '_recog_erp_ar_diff_waves_lpfilt_roi.erp' ...
        };
    roi_erpname_suffixes = { ...
        '_recog_erp_ar_diff_waves_roi', ...
        '_recog_erp_ar_diff_waves_lpfilt_roi' ...
        };
    roi_list_paths = { ...
        fullfile(output_dir, 'Measurement_ERP_List_recog_roi.txt'), ...
        fullfile(output_dir, 'Measurement_ERP_List_lpfilt_recog_roi.txt') ...
        };

    for li = 1:numel(roi_input_suffixes)
        fid = fopen(roi_list_paths{li}, 'w');
        for i = 1:length(SUB)
            Subject_Path = fullfile(data_dir, SUB{i});
            inName = [SUB{i} roi_input_suffixes{li}];
            inFile = fullfile(Subject_Path, inName);
            if exist(inFile, 'file') ~= 2
                error('ERP file not found: %s', inFile);
            end

            ERP = pop_loaderp('filename', inName, 'filepath', Subject_Path);
            base_nchan = ERP.nchan;
            if base_nchan ~= numel(labels_all)
                error('Channel count mismatch for %s: expected %d, got %d', ...
                    SUB{i}, numel(labels_all), base_nchan);
            end
            if max([chan_fn400 chan_lpc]) > base_nchan
                error('Channel index exceeds channel count for %s', SUB{i});
            end

            ERP.bindata(base_nchan+1,:,:) = mean(ERP.bindata(chan_fn400,:,:), 1, 'omitnan');
            ERP.bindata(base_nchan+2,:,:) = mean(ERP.bindata(chan_lpc,:,:), 1, 'omitnan');

            if isfield(ERP, 'binerror') && size(ERP.binerror,1) == base_nchan
                ERP.binerror(base_nchan+1,:,:) = mean(ERP.binerror(chan_fn400,:,:), 1, 'omitnan');
                ERP.binerror(base_nchan+2,:,:) = mean(ERP.binerror(chan_lpc,:,:), 1, 'omitnan');
            end
            if isfield(ERP, 'binsem') && size(ERP.binsem,1) == base_nchan
                ERP.binsem(base_nchan+1,:,:) = mean(ERP.binsem(chan_fn400,:,:), 1, 'omitnan');
                ERP.binsem(base_nchan+2,:,:) = mean(ERP.binsem(chan_lpc,:,:), 1, 'omitnan');
            end

            newLoc = ERP.chanlocs(1);
            if isfield(newLoc, 'X'); newLoc.X = NaN; end
            if isfield(newLoc, 'Y'); newLoc.Y = NaN; end
            if isfield(newLoc, 'Z'); newLoc.Z = NaN; end
            if isfield(newLoc, 'theta'); newLoc.theta = []; end
            if isfield(newLoc, 'radius'); newLoc.radius = []; end
            if isfield(newLoc, 'sph_theta'); newLoc.sph_theta = []; end
            if isfield(newLoc, 'sph_phi'); newLoc.sph_phi = []; end
            if isfield(newLoc, 'sph_radius'); newLoc.sph_radius = []; end
            if isfield(newLoc, 'type'); newLoc.type = ''; end
            if isfield(newLoc, 'urchan'); newLoc.urchan = base_nchan + 1; end
            newLoc.labels = roi_label_fn400;
            ERP.chanlocs(base_nchan+1) = newLoc;

            if isfield(newLoc, 'urchan'); newLoc.urchan = base_nchan + 2; end
            newLoc.labels = roi_label_lpc;
            ERP.chanlocs(base_nchan+2) = newLoc;

            ERP.nchan = base_nchan + 2;

            outFile = fullfile(roi_output_dir, [SUB{i} roi_output_suffixes{li}]);
            ERP = pop_savemyerp(ERP, 'erpname', [SUB{i} roi_erpname_suffixes{li}], 'filename', outFile);
            fprintf(fid, '%s\n', outFile);
        end
        fclose(fid);
    end
end

%*************************************************************************************************************************************

if do_roi_average
    ERPset_list_unfilt_roi = roi_list_paths{1};
    ERPset_list_lpfilt_roi = roi_list_paths{2};
end

%*************************************************************************************************************************************

%Difference waveform measurements on averaged ERP waveforms without a low-pass filter applied
%Create a text file containing a list of unfiltered ERPsets and their file locations to measure mean amplitude and 50% area latency from
ERPset_list_unfilt = fullfile(output_dir, 'Measurement_ERP_List_recog.txt');
fid = fopen(ERPset_list_unfilt, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_recog_erp_ar_diff_waves.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Measure mean amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_fn400, diffbin, chan_fn400, 'Baseline', baselinecorr, 'Measure', 'meanbl',... 
    'Filename', fullfile(output_dir, 'Mean_Amplitude_Diff_Waves_recog.txt'), 'Binlabel', 'on', 'FileFormat', 'wide',... 
    'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_fn400, diffbin, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure', 'meanbl',... 
        'Filename', fullfile(output_dir, 'Mean_Amplitude_Diff_Waves_recog_FN400_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide',... 
        'InterpFactor',  1,  'Resolution', 3);
end

%Measure 50% area latency (positive area only) using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_fn400, diffbin, chan_fn400, 'Baseline', baselinecorr, 'Measure', 'fareaplat',... 
    'Afraction', 0.5, 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename', fullfile(output_dir, '50%_Area_Latency_Diff_Waves_recog.txt'),...
    'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_fn400, diffbin, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure', 'fareaplat',... 
        'Afraction', 0.5, 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename', fullfile(output_dir, '50%_Area_Latency_Diff_Waves_recog_FN400_ROI.txt'),...
        'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%*************************************************************************************************************************************

%Difference waveform measurements on averaged ERP waveforms with a low-pass filter applied

%Create a text file containing a list of low-pass filtered ERPsets and their file locations to measure peak amplitude, peak latency and onset latency (50% peak latency) from
ERPset_list_lpfilt = fullfile(output_dir, 'Measurement_ERP_List_lpfilt_recog.txt');
fid = fopen(ERPset_list_lpfilt, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_recog_erp_ar_diff_waves_lpfilt.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Measure (local) peak amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_fn400, diffbin, chan_fn400, 'Baseline', baselinecorr, 'Measure','peakampbl',... 
    'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',...
    fullfile(output_dir, 'Peak_Amplitude_Diff_Waves_recog.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_fn400, diffbin, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure','peakampbl',... 
        'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',...
        fullfile(output_dir, 'Peak_Amplitude_Diff_Waves_recog_FN400_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%Measure (local) peak latency using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_fn400, diffbin, chan_fn400, 'Baseline', baselinecorr, 'Measure','peaklatbl',... 
    'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',... 
    fullfile(output_dir, 'Peak_Latency_Diff_Waves_recog.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_fn400, diffbin, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure','peaklatbl',... 
        'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',... 
        fullfile(output_dir, 'Peak_Latency_Diff_Waves_recog_FN400_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%Measure onset latency (50% peak latency) using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_onsetlat_fn400, diffbin, chan_fn400, 'Baseline', baselinecorr, 'Measure','fpeaklat', 'Neighborhood',  3,... 
    'Peakpolarity', 'positive', 'Afraction', 0.5, 'Neighborhood',  3, 'Peakreplace', 'absolute', 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename',...
    fullfile(output_dir, 'Onset_Latency_Diff_Waves_recog.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_onsetlat_fn400, diffbin, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure','fpeaklat', 'Neighborhood',  3,... 
        'Peakpolarity', 'positive', 'Afraction', 0.5, 'Neighborhood',  3, 'Peakreplace', 'absolute', 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename',...
        fullfile(output_dir, 'Onset_Latency_Diff_Waves_recog_FN400_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%*************************************************************************************************************************************

%Parent waveform measurements on averaged ERP waveforms without a low-pass filter applied

%Measure mean amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_fn400, parentbins, chan_fn400, 'Baseline', baselinecorr, 'Measure', 'meanbl', 'Filename',... 
    fullfile(output_dir, 'Mean_Amplitude_Parent_Waves_recog.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_fn400, parentbins, roi_chan_fn400, 'Baseline', baselinecorr, 'Measure', 'meanbl', 'Filename',... 
        fullfile(output_dir, 'Mean_Amplitude_Parent_Waves_recog_FN400_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%*************************************************************************************************************************************

%% LPC measurements

%Difference waveform measurements on averaged ERP waveforms without a low-pass filter applied
ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_lpc, diffbin, chan_lpc, 'Baseline', baselinecorr, 'Measure', 'meanbl',... 
    'Filename', fullfile(output_dir, 'Mean_Amplitude_Diff_Waves_recog_LPC.txt'), 'Binlabel', 'on', 'FileFormat', 'wide',... 
    'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_lpc, diffbin, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure', 'meanbl',... 
        'Filename', fullfile(output_dir, 'Mean_Amplitude_Diff_Waves_recog_LPC_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide',... 
        'InterpFactor',  1,  'Resolution', 3);
end

ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_lpc, diffbin, chan_lpc, 'Baseline', baselinecorr, 'Measure', 'fareaplat',... 
    'Afraction', 0.5, 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename', fullfile(output_dir, '50%_Area_Latency_Diff_Waves_recog_LPC.txt'),...
    'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_lpc, diffbin, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure', 'fareaplat',... 
        'Afraction', 0.5, 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename', fullfile(output_dir, '50%_Area_Latency_Diff_Waves_recog_LPC_ROI.txt'),...
        'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%Difference waveform measurements on averaged ERP waveforms with a low-pass filter applied
ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_lpc, diffbin, chan_lpc, 'Baseline', baselinecorr, 'Measure','peakampbl',... 
    'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',...
    fullfile(output_dir, 'Peak_Amplitude_Diff_Waves_recog_LPC.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_lpc, diffbin, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure','peakampbl',... 
        'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',...
        fullfile(output_dir, 'Peak_Amplitude_Diff_Waves_recog_LPC_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_lpc, diffbin, chan_lpc, 'Baseline', baselinecorr, 'Measure','peaklatbl',... 
    'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',... 
    fullfile(output_dir, 'Peak_Latency_Diff_Waves_recog_LPC.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_lpc, diffbin, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure','peaklatbl',... 
        'Peakpolarity', 'positive', 'Neighborhood',  3, 'PeakReplace', 'absolute', 'Filename',... 
        fullfile(output_dir, 'Peak_Latency_Diff_Waves_recog_LPC_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

ALLERP = pop_geterpvalues( ERPset_list_lpfilt, timewindow_onsetlat_lpc, diffbin, chan_lpc, 'Baseline', baselinecorr, 'Measure','fpeaklat', 'Neighborhood',  3,... 
    'Peakpolarity', 'positive', 'Afraction', 0.5, 'Neighborhood',  3, 'Peakreplace', 'absolute', 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename',...
    fullfile(output_dir, 'Onset_Latency_Diff_Waves_recog_LPC.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_lpfilt_roi, timewindow_onsetlat_lpc, diffbin, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure','fpeaklat', 'Neighborhood',  3,... 
        'Peakpolarity', 'positive', 'Afraction', 0.5, 'Neighborhood',  3, 'Peakreplace', 'absolute', 'PeakOnset',  1, 'Fracreplace', 'NaN', 'Filename',...
        fullfile(output_dir, 'Onset_Latency_Diff_Waves_recog_LPC_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%Parent waveform measurements on averaged ERP waveforms without a low-pass filter applied
ALLERP = pop_geterpvalues( ERPset_list_unfilt, timewindow_lpc, parentbins, chan_lpc, 'Baseline', baselinecorr, 'Measure', 'meanbl', 'Filename',... 
    fullfile(output_dir, 'Mean_Amplitude_Parent_Waves_recog_LPC.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
if do_roi_average
    ALLERP = pop_geterpvalues( ERPset_list_unfilt_roi, timewindow_lpc, parentbins, roi_chan_lpc, 'Baseline', baselinecorr, 'Measure', 'meanbl', 'Filename',... 
        fullfile(output_dir, 'Mean_Amplitude_Parent_Waves_recog_LPC_ROI.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);
end

%*************************************************************************************************************************************
