%Script #12
%Operates on individual subject data
%Uses the output from Script #7: Average_ERPs.m
%This script uses the individual subject averaged ERP waveforms from Script #7, and measures the mean amplitude, peak amplitude, peak latency, 50% area latency, and onset latency (50% peak latency) 
%during the time window of the component, and saves a separate text file for each measurement in the ERP Measurements folder. 
%Note that based on their respective susceptibility to high frequency noise, some measurements (e.g., mean amplitude, 50% area latency) are calculated on the averaged ERP waveforms without a low-pass filter applied, 
%whereas other measurements (e.g., peak amplitude, peak latency, and onset latency) are calculated from waveforms that have been low-pass filtered. 

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject ERP files are stored in the sibling data directory (eeg/audi/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'outputs');
erp_output_dir = fullfile(output_dir, 'ERP_Measurements');
data_dir = fileparts(project_dir);

if ~exist(erp_output_dir, 'dir')
    mkdir(erp_output_dir);
end

DIR = data_dir;

fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP12 Measure ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
%sub-010质量不好
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', 'sub-008', 'sub-009'};   
task ='audi';

%*************************************************************************************************************************************

%Set measurement time window for measuring mean amplitude, peak amplitude, peak latency, and 50% area latency in milliseconds (e.g., 300 to 600 ms)
timewindow = [200 500];

%Set measurement time window for measuring onset latency (50% peak latency) in milliseconds (e.g., 200 to 600 ms)
timewindow_onsetlat = [100 500];

%Set EEG channel(s) to measure the components  (FCz= 26 Cz=35, CPz=44, Pz=53)
%画Cz=35
chan = [35];

%Set difference wave bin(s) for measurement
diffbin = [3];  

%Set parent wave bins for measurement
parentbins = [1 2];  

%Set baseline correction period for measurement
baselinecorr = [-200 0]; 

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%*************************************************************************************************************************************

%Difference waveform measurements on averaged ERP waveforms without a low-pass filter applied

%Create a text file containing a list of unfiltered ERPsets and their file locations to measure mean amplitude and 50% area latency from
ERPset_list = fullfile(erp_output_dir, 'Measurement_ERP_List_audi.txt');
fid = fopen(ERPset_list, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Measure mean amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow, diffbin, chan, 'Baseline', baselinecorr, 'Measure', 'meanbl',... 
    'Filename', fullfile(erp_output_dir, 'Mean_Amplitude_Diff_Waves_audi.txt'), 'Binlabel', 'on', 'FileFormat', 'wide',... 
    'InterpFactor',  1,  'Resolution', 3);

%Measure 50% area latency (positive area only) using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow, diffbin, chan, 'Baseline', baselinecorr, 'Measure', 'fareaplat',... 
    'Afraction', 0.5, 'PeakOnset', 1, 'Fracreplace', 'NaN', 'Filename', fullfile(erp_output_dir, '50%_Area_Latency_Diff_Waves_audi.txt'),...
    'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor',  1,  'Resolution', 3);

%*************************************************************************************************************************************

%Difference waveform measurements on averaged ERP waveforms with a low-pass filter applied

%Create a text file containing a list of low-pass filtered ERPsets and their file locations to measure peak amplitude, peak latency and onset latency (50% peak latency) from
ERPset_list = fullfile(erp_output_dir, 'Measurement_ERP_List_lpfilt_audi.txt');
fid = fopen(ERPset_list, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves_lpfilt.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Measure (local) peak amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow, diffbin, chan, 'Baseline', baselinecorr, 'Measure', 'peakampbl',... 
    'Peakpolarity', 'positive', 'Neighborhood', 3, 'PeakReplace', 'absolute', 'Filename', ...
    fullfile(erp_output_dir, 'Peak_Amplitude_Diff_Waves_audi.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor', 1, 'Resolution', 3);

%Measure (local) peak latency using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow, diffbin, chan, 'Baseline', baselinecorr, 'Measure', 'peaklatbl',... 
    'Peakpolarity', 'positive', 'Neighborhood', 3, 'PeakReplace', 'absolute', 'Filename', ... 
    fullfile(erp_output_dir, 'Peak_Latency_Diff_Waves_audi.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor', 1, 'Resolution', 3);

%Measure onset latency (50% peak latency) using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow_onsetlat, diffbin, chan, 'Baseline', baselinecorr, 'Measure', 'fpeaklat', 'Neighborhood', 3,... 
    'Peakpolarity', 'positive', 'Afraction', 0.5, 'Neighborhood', 3, 'Peakreplace', 'absolute', 'PeakOnset', 1, 'Fracreplace', 'NaN', 'Filename', ...
    fullfile(erp_output_dir, 'Onset_Latency_Diff_Waves_audi.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor', 1, 'Resolution', 3);

%*************************************************************************************************************************************

%Parent waveform measurements on averaged ERP waveforms without a low-pass filter applied

ERPset_list = fullfile(erp_output_dir, 'Measurement_ERP_List_audi.txt');

%Measure mean amplitude using the time window, channel(s), and bin(s) specified above
ALLERP = pop_geterpvalues(ERPset_list, timewindow, parentbins, chan, 'Baseline', baselinecorr, 'Measure', 'meanbl', 'Filename',... 
    fullfile(erp_output_dir, 'Mean_Amplitude_Parent_Waves_audi.txt'), 'Binlabel', 'on', 'FileFormat', 'wide', 'InterpFactor', 1, 'Resolution', 3);

%*************************************************************************************************************************************
