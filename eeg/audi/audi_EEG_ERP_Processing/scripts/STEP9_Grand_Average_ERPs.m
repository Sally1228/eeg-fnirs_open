%Script #9
%Operates on individual subject data
%Uses the output from Script #7: Average_ERPs.m
%This script uses the individual subject averaged ERP waveforms from Script #7 to create grand average ERP waveforms across participants both with and without a low-pass filter applied. 

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject ERP files are stored in the sibling data directory (eeg/audi/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'outputs');
ga_output_dir = fullfile(output_dir, 'Grand_Average_ERPs');
data_dir = fileparts(project_dir);

if ~exist(ga_output_dir, 'dir')
    mkdir(ga_output_dir);
end


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP9 Grand average ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', 'sub-008', 'sub-009', 'sub-010'};   
task ='audi';


%*************************************************************************************************************************************

%Create grand average ERP waveforms from individual subject ERPs without low-pass filter applied 

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Create a text file containing a list of ERPsets and their file locations to include in the grand average ERP waveforms
ERPset_list = fullfile(ga_output_dir, 'GA_audi_erp_ar_diff_waves.txt');
fid = fopen(ERPset_list, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Create a grand average ERP waveform
ERP = pop_gaverager( ERPset_list , 'ExcludeNullBin', 'on', 'SEM', 'on' );
ERP = pop_savemyerp(ERP, 'erpname', 'GA_audi_erp_ar_diff_waves', 'filename', 'GA_audi_erp_ar_diff_waves.erp', 'filepath', ga_output_dir, 'Warning', 'off');

%*************************************************************************************************************************************

%Create grand average ERP waveforms from individual subject ERPs with a low-pass filter applied

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Create a text file containing a list of low-pass filtered ERPsets and their file locations to include in the grand average ERP waveforms
ERPset_list = fullfile(ga_output_dir, 'GA_audi_erp_ar_diff_waves_lpfilt.txt');
fid = fopen(ERPset_list, 'w');
    for i = 1:length(SUB)
        Subject_Path = fullfile(data_dir, SUB{i});
        erppath = fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves_lpfilt.erp']);
        fprintf(fid,'%s\n', erppath);
    end
fclose(fid);

%Create a grand average ERP waveform
ERP = pop_gaverager( ERPset_list , 'ExcludeNullBin', 'on', 'SEM', 'on' );
ERP = pop_savemyerp(ERP, 'erpname', 'GA_audi_erp_ar_diff_waves_lpfilt', 'filename', 'GA_audi_erp_ar_diff_waves_lpfilt.erp', 'filepath', ga_output_dir, 'Warning', 'off');

%*************************************************************************************************************************************
