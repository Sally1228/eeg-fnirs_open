%Script #5
%Operates on individual subject data
%Uses the output from Script #4: Remove_ICA_Components.m
%This script loads the semi-continuous ICA-corrected EEG data file from Script #4, creates an Event List containing a record of all event codes and their timing, assigns events to bins using Binlister, epochs the EEG, and performs baseline correction.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/cate/sub-XXX).
% script_dir =  '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/cate/cate_EEG_ERP_Processing/scripts'
script_dir = fileparts(mfilename('fullpath'));

project_dir= fileparts(script_dir);
input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP5 Event list, binning, epoch\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-001','sub-002','sub-003', 'sub-004', 'sub-005', 'sub-008', 'sub-009' };
task ='cate';

%**********************************************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the semi-continuous ICA-corrected EEG data file outputted from Script #4 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr'], 'overwrite','on','gui', 'off'); 

    % 创建事件列表 Create EEG Event List containing a record of all event codes and their timing
    EEG  = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, ...
        'Eventlist', fullfile(Subject_Path, [SUB{i} '_cate_Eventlist.txt'])); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist.set']), 'overwrite', 'on', 'gui', 'off');

    %Assign events to bins with Binlister; an individual trial may be assigned to more than one bin (bin assignments can be reviewed in each subject's cate_Eventlist_Bins.txt file)
    EEG  = pop_binlister(EEG, 'BDF', fullfile(input_dir, 'BDF_cate.txt'), ...
        'ExportEL', fullfile(Subject_Path, [SUB{i} '_cate_Eventlist_Bins.txt']), 'IndexEL', 1, 'SendEL2', 'EEG&Text', 'UpdateEEG', 'on', 'Voutput', 'EEG');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins.set']), 'overwrite', 'on', 'gui', 'off'); 

    %Epoch the EEG 60s任务，休息15s
    %基线0-2s，任务0-60s，先不做基线矫正
    EEG = pop_epochbin(EEG, [-2000 60000], 'none');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch.set']), 'overwrite', 'on', 'gui', 'off'); 
    close all;
    
%End subject loop
end

%**********************************************************************************************************************************************************************
