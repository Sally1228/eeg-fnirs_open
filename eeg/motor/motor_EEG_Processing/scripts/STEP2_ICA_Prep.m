%Script #2
%Operates on individual subject data
%Uses the output from Script #1: Import_Raw_EEG_Shift_Reref_DS_Hpfilt.m
%This script loads the outputted continuous EEG data file from Script #1, removes segments of EEG during the break periods in between trial blocks, and
%removes especially noisy segments of EEG during the trial blocks to prepare the data for ICA. Note that the goal of this stage of processing is to remove 
%particularly noisy segments of data; a more thorough rejection of artifacts will be performed later on the epoched data.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/motor/sub-XXX).
% script_dir =  '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/motor/motor_EEG_ERP_Processing/scripts'
script_dir = fileparts(mfilename('fullpath'));

project_dir= fileparts(script_dir);
input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);

fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP2 ICA prep\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-001','sub-002','sub-003', 'sub-004', 'sub-005', 'sub-006','sub-007', 'sub-008', 'sub-009', 'sub-010' };
task ='motor';

%*************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the continuous EEG data file outputted from Script #1 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_notch.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_shifted_ds_reref_ucbip_hpfilt'], 'gui', 'off'); 

    %Load parameters for rejecting especially noisy segments of EEG during trial blocks from Excel file ICA_Prep_Values_P3.xls. Default parameters can be used initially but may need 
    % to be modified for a given participant on the basis of visual inspection of the data.
    [ndata, text, alldata] = xlsread(fullfile(input_dir, 'ICA_Prep_Values_motor.xlsx')); 
        for j = 1:length(alldata)           
            if isequal(SUB{i},num2str(alldata{j,1}));
                AmpthValue = alldata{j,2};
                WindowValue = alldata{j,3};
                StepValue = alldata{j,4};
            end
        end

    %Delete segments of the EEG exceeding the thresholds defined above
    %rewiew=on 可以用于检查
    %滑动窗口越短，对瞬时尖峰更敏感；越长，更偏向抓持续性大幅度漂移或肌电。
    %窗口滑动步长ms 越小扫描越细但更慢，也更容易把一段伪迹标得更连续。
    EEG = pop_continuousartdet( EEG, 'ampth', AmpthValue, 'winms', WindowValue, 'stepms', StepValue, 'chanArray', 1:62, 'review', 'off');        
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_shifted_ds_reref_ucbip_hpfilt_ica_pre2'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_pre2.set']), 'overwrite', 'on', 'gui', 'off'); 

%End subject loop
end

%*************************************************************************************************************************************
