%Script #4
%Operates on individual subject data
%Uses the output from Script #3: Run_ICA.m
%This script loads the outputted semi-continuous EEG data file containing the ICA weights from Script #3, loads the list of ICA component(s) from the ICA_Components_P3.xlsx Excel file, and removes the component(s) from the EEG.
%Note that if ICA weights were re-computed on the data, the component(s) to remove will need to be updated in the Excel file to match the new components (see Script #3: Run_ICA.m for further details).

%% Config
close all; clearvars;

% script_dir =  '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/motor/motor_EEG_ERP_Processing/scripts'
script_dir = fileparts(mfilename('fullpath'));

project_dir= fileparts(script_dir);
input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);

fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP4 Remove ICA components\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-001','sub-002','sub-003', 'sub-004', 'sub-005', 'sub-006','sub-007', 'sub-008', 'sub-009', 'sub-010' };
task ='motor';

%**********************************************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the continuous EEG data file containing the ICA weights outputted from Script #3 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted'],'overwrite','on', 'gui','off'); 

    %Load list of ICA component(s) corresponding to ocular artifacts from Excel file ICA_Components_P3.xlsx
    [ndata, text, alldata] = xlsread(fullfile(input_dir, 'ICA_Components_motor.xlsx')); 
    MaxNumComponents = size(alldata, 2);
        for j = 1:length(alldata)
            if isequal(SUB{i}, num2str(alldata{j,1}));
                NumComponents = 0;
                for k = 2:MaxNumComponents
                    if ~isnan(alldata{j,k});
                        NumComponents = NumComponents+1;
                    end
                    Components = [alldata{j,(2:(NumComponents+1))}];
                end
            end
        end

    %Perform ocular correction by removing the ICA component(s) specified above
    EEG = pop_subcomp( EEG, [Components], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr.set']), 'overwrite', 'on', 'gui', 'off'); 
    
%End subject loop
end

%**********************************************************************************************************************************************************************
