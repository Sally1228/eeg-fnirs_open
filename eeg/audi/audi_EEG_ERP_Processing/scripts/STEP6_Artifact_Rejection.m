%Script #6 
% 注意:eeg被标记，但是没有被剔除，在后续 pop_averager 时设置"排除被标记的试次"
% 需要文档：Interpolate_Channels_audi、AR_Parameters_for_SVT_CRAP_audi、AR_Parameters_for_MW_CRAP_audi
%Operates on individual subject data
%Uses the output from Script #5: Elist_Bin_Epoch.m
%This script loads the epoched EEG data file from Script #5, interpolates bad channels listed in Excel file Interpolate_Channels_P3.xls, and performs artifact rejection to remove noisy segments of EEG, 
%segments containing eyeblinks or eye movements during the time of the stimulus (i.e., resulting in a change in sensory input on that trial), and segments containing uncorrected residual eye movements throughout the epoch
%using the parameters tailored to an individual subject's data listed in the corresponding Excel file for that artifact.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/audi/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP6 Artifact rejection\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', 'sub-008', 'sub-009', 'sub-010'};   
task ='audi';

%**********************************************************************************************************************************************************************
%预先设置的参数
%Load the Excel file with the list of channels to interpolate for each subject 
[ndata1, text1, alldata1] = xlsread(fullfile(input_dir, 'Interpolate_Channels_audi.xlsx'));
%Load the Excel file with the list of thresholds and parameters for identifying C.R.A.P. with the simple voltage threshold algorithm for each subject 
[ndata2, text2, alldata2] = xlsread(fullfile(input_dir, 'AR_Parameters_for_SVT_CRAP_audi.xlsx'));
%Load the Excel file with the list of thresholds and parameters for identifying C.R.A.P. with the moving window peak-to-peak algorithm for each subject 
[ndata3, text3, alldata3] = xlsread(fullfile(input_dir, 'AR_Parameters_for_MW_CRAP_audi.xlsx'));

%*************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the epoched EEG data file outputted from Script #5 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch'], 'overwrite','on', 'gui', 'off'); 

    %Interpolate channel(s) specified in Excel file Interpolate_Channels_audi.xls; any channels without channel locations (e.g., the eye channels) should not be included in the interpolation process and are listed in ignored channels
    %EEG channels that will later be used for measurement of the ERPs should not be interpolated

    %补充EEG.chanlocs.urchan，否认下面会报错
    if ~isfield(EEG.chanlocs, 'urchan')
        for m = 1:EEG.nbchan
            EEG.chanlocs(m).urchan = m;
        end
    end
    
    %耳电极不参与插补
    ignored_channels = [61 62];        
    DimensionsOfFile1 = size(alldata1);
    for j = 1:DimensionsOfFile1(1);
        if isequal(SUB{i},num2str(alldata1{j,1}));
           badchans = (alldata1{j,2});
           if ~isequal(badchans,'none') | ~isempty(badchans)
           	  if ~isnumeric(badchans)
                 badchans = str2num(badchans);
              end
              EEG  = pop_erplabInterpolateElectrodes( EEG , 'displayEEG',  0, 'ignoreChannels',  ignored_channels, 'interpolationMethod', 'spherical', 'replaceChannels', badchans);
           end
           [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp'], ...
               'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp.set']), 'overwrite', 'on', 'gui', 'off'); 
        end
    end

    %Identify segments of EEG with C.R.A.P. artifacts using the simple voltage threshold algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile2 = size(alldata2);
    for j = 1:DimensionsOfFile2(1)
        if isequal(SUB{i},num2str(alldata2{j,1}));
            if isequal(alldata2{j,2}, 'default')
                Channels = 1:62;
            else
                Channels = str2num(alldata2{j,2});
            end
            ThresholdMinimum = alldata2{j,3};
            ThresholdMaximum = alldata2{j,4};
            TimeWindowMinimum = alldata2{j,5};
            TimeWindowMaximum = alldata2{j,6};
        end
    end

    EEG  = pop_artextval( EEG , 'Channel',  Channels, 'Flag', [1 2], 'Threshold', [ThresholdMinimum ThresholdMaximum], 'Twindow', [TimeWindowMinimum  TimeWindowMaximum] ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_SVT'], 'overwrite','on','gui', 'off'); 

    %Identify segments of EEG with C.R.A.P. artifacts using the moving window peak-to-peak algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile3 = size(alldata3);
    for j = 1:DimensionsOfFile3(1)
        if isequal(SUB{i},num2str(alldata3{j,1}));
            if isequal(alldata3{j,2}, 'default')
                Channels = 6:60; %剔除耳电极、前额
            else
                Channels = str2num(alldata3{j,2});
            end
            Threshold = alldata3{j,3};
            TimeWindowMinimum = alldata3{j,4};
            TimeWindowMaximum = alldata3{j,5};
            WindowSize = alldata3{j,6};
            WindowStep = alldata3{j,7};
        end
    end

    EEG  = pop_artmwppth( EEG , 'Channel',  Channels, 'Flag', [1 3], 'Threshold', Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize', WindowSize, 'Windowstep', WindowStep ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_SVT_MW1'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set']), 'overwrite', 'on', 'gui', 'off'); 

end

%*************************************************************************************************************************************
