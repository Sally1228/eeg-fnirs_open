%Script #1
%Operates on individual subject data
%This script loads the raw continuous EEG data in .set EEGLAB file format, shifts the stimulus event codes forward in time to account for the LCD monitor 
%delay (26 ms on our monitor, as measured with a photosensor), downsamples the data to 256 Hz to speed data processing time, references to the average of 
%P9 and P10, creates bipolar HEOG and VEOG channels, adds channel location information, removes the DC offsets, and applies a high-pass filter. 


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
fprintf('[%s] STEP1 Import raw EEG\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-001','sub-002','sub-003', 'sub-004', 'sub-005', 'sub-007', 'sub-008', 'sub-009', 'sub-010' };
task ='cate';


%***********************************************************************************************************************************************
%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the raw continuous EEG data file in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_task-' task '_eeg'], 'overwrite', 'on', 'gui', 'off'); 

    %Downsample from the recorded sampling rate of 1000 Hz to 256 Hz to speed data processing (automatically applies the appropriate low-pass anti-aliasing filter)
    EEG = pop_resample( EEG, 256);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_task-' task '_eeg_ds'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds.set']), 'overwrite', 'on', 'gui', 'off'); 

    %Rereference to the average of A1 and A2 以耳电极平均作为参考/平均参考
    EEG = pop_eegchanoperator(EEG, fullfile(input_dir, 'Rereference_Add_Uncorrected_A1A2.txt'));
    %EEG = pop_eegchanoperator(EEG, fullfile(Current_File_Path, 'Rereference_Add_Uncorrected_av_cate.txt'));

    %Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
    EEG = pop_chanedit(EEG, 'lookup', fullfile(input_dir, 'standard-10-5-cap385.elp'));
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref.set']), 'overwrite', 'on', 'gui', 'off');

    %Remove DC offsets and apply a high-pass filter (non-causal Butterworth impulse response function, 0.1 Hz half-amplitude cut-off, 12 dB/oct roll-off)
    EEG  = pop_basicfilter( EEG,  1:62 , 'Boundary', 'boundary', 'Cutoff',  0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' );
    %50Hz工频滤波 luck推荐用cleanline，发现滤不干净。先用常规方法。
    %EEG = pop_cleanline(EEG, 'SignalType','Channels','ChanCompIndices',1:EEG.nbchan, 'LineFrequencies',[50 100], 'ScanForLines',true, 'LineAlpha',0.01, 'Bandwidth',2,'SlidingWinLength',4, 'SlidingWinStep',2, 'SmoothingFactor',100, 'ComputeSpectralPower',true,'PlotFigures',false,'VerbosityLevel',1);
    EEG = pop_eegfiltnew(EEG, 49, 51, [], 1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_notch'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_notch.set']), 'overwrite', 'on', 'gui', 'off');

%End subject loop
end

%***********************************************************************************************************************************************
