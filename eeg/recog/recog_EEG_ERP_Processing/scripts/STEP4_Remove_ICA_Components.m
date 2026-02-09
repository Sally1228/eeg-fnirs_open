%Script #4
%Operates on individual subject data
%Uses the output from Script #3: Run_ICA.m
%This script loads the outputted semi-continuous EEG data file containing the ICA weights from Script #3, loads the list of ICA component(s) from the ICA_Components_P3.xlsx Excel file, and removes the component(s) from the EEG.
%Note that if ICA weights were re-computed on the data, the component(s) to remove will need to be updated in the Excel file to match the new components (see Script #3: Run_ICA.m for further details).

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/recog/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP4 Remove ICA components\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-008', 'sub-009', 'sub-010'};   
task ='recog';

%**********************************************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %Define subject path based on study directory and subject ID of current subject
     Subject_Path = fullfile(data_dir, SUB{i});

    %Load the continuous EEG data file containing the ICA weights outputted from Script #3 in .set EEGLAB file format
    EEG = pop_loadset('filename',[SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted.set'],'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted'],'overwrite','on', 'gui','off'); 

    %Load list of ICA component(s) corresponding to ocular artifacts from Excel file ICA_Components_P3.xlsx
    [ndata, text, alldata] = xlsread(fullfile(input_dir, 'ICA_Components_recog.xlsx')); 
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
    %剔除眼动和工频50Hz
    EEG = pop_subcomp( EEG, [Components], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr.set']), 'overwrite','on', 'gui','off'); 
    
    %Create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2) from the ICA corrected data; the original uncorrected HEOG and VEOG channels are retained for later artifact detection procedures
    %EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'Add_Corrected_Bipolars_P3.txt']);    
    %源代码报错，修改成递归方式
    %EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'Add_Corrected_Bipolars_P3_recursive.txt']);
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_P3_shifted_ds_reref_ucbip_hpfilt_ica_corr_35chan'],'gui', 'off'); 
    
    %Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
    %EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'standard-10-5-cap385.elp']);
    %cbip = corrected bipolar
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_P3_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip'], 'savenew', [Subject_Path SUB{i} '_P3_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip.set'], 'gui', 'off'); 

%End subject loop
end

%**********************************************************************************************************************************************************************
