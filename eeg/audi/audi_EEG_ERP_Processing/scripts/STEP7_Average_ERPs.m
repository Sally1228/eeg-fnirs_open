%Script #7
%Operates on individual subject data
%Uses the output from Script #6: Artifact_Rejection.m
%This script loads the epoched and artifact rejected EEG data from Script #6, creates an averaged ERP waveform, calculates the percentage of trials rejected for artifacts (in total and per bin) 
%and saves the information to a .csv file in each subject's data folder, calculates ERP difference waveforms between conditions, and creates low-pass filtered versions of the ERP waveforms.


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
fprintf('[%s] STEP7 Average ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', 'sub-008', 'sub-009', 'sub-010'};   
task ='audi';


%**********************************************************************************************************************************************************************

%Create averaged ERP waveforms

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));
  
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the epoched and artifact rejected EEG data file outputted from Script #6 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set'], 'filepath', Subject_Path);

    %Create an averaged ERP waveform
    %'Criterion','good'：只平均 good trials，也就是那些没被你前面伪迹检测标记为 artifact、未被剔除的 trial。
    % 'ExcludeBoundary','on'：如果某个 epoch 跨越了 boundary 事件（比如删段、拼接造成的断点），就不把它纳入平均，避免边界伪迹污染平均波形。
    % 'SEM','on'：同时计算 SEM（standard error of the mean，均值标准误），方便后面画误差带或报告变异度。
    ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on');
    ERP = pop_savemyerp(ERP, 'erpname', [SUB{i} '_audi_erp_ar'], 'filename', fullfile(Subject_Path, [SUB{i} '_audi_erp_ar.erp']));
    
    %Apply a low-pass filter (non-causal Butterworth impulse response function, 20 Hz half-amplitude cut-off, 48 dB/oct roll-off) to the ERP waveforms
    %低通20Hz
    ERP = pop_filterp( ERP,  1:62 , 'Cutoff',  20, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  8 );
    ERP = pop_savemyerp(ERP, 'erpname', [SUB{i} '_audi_erp_ar_lpfilt'], 'filename', fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_lpfilt.erp']));

    %Calculate the percentage of trials that were rejected in each bin 
    accepted = ERP.ntrials.accepted;
    rejected= ERP.ntrials.rejected;
    percent_rejected= rejected./(accepted + rejected)*100;
    
    %Calculate the total percentage of trials rejected across all trial types (first two bins)
    total_accepted = accepted(1) + accepted(2);
    total_rejected= rejected(1)+ rejected(2);
    total_percent_rejected= total_rejected./(total_accepted + total_rejected)*100; 
    
    %Save the percentage of trials rejected (in total and per bin) to a .csv file 
    fid = fopen(fullfile(data_dir, SUB{i}, [SUB{i} '_AR_Percentages_audi.csv']), 'w');
    fprintf(fid, 'SubID,Bin,Accepted,Rejected,Total Percent Rejected\n');
    fprintf(fid, '%s,%s,%d,%d,%.2f\n', SUB{i}, 'Total', total_accepted, total_rejected, total_percent_rejected);
    bins = strrep(ERP.bindescr,', ',' - ');
    for b = 1:length(bins)
        fprintf(fid, ',%s,%d,%d,%.2f\n', bins{b}, accepted(b), rejected(b), percent_rejected(b));
    end
    fclose(fid);
    
%End subject loop
end

%**********************************************************************************************************************************************************************

%Create difference waveforms 

%Loop through each subject listed in SUB
for i = 1:length(SUB)
   
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(DIR, SUB{i]);        
    %Load averaged ERP waveform (without the 20 Hz low-pass filter) 
    ERP = pop_loaderp('filename', [SUB{i} '_audi_erp_ar_lpfilt.erp'], 'filepath', Subject_Path);                                                                                                                                                                                                                                                                                                                                       

    %Create ERP difference waveforms between conditions
    ERP = pop_binoperator(ERP, fullfile(input_dir, 'audi_Diff_Wave.txt'));
    ERP = pop_savemyerp(ERP, 'erpname', [SUB{i} '_audi_erp_ar_diff_waves'], 'filename', fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves.erp']));

    %Apply a low-pass filter (non-causal Butterworth impulse response function, 20 Hz half-amplitude cut-off, 48 dB/oct roll-off) to the difference waveforms
    ERP = pop_filterp( ERP,  1:62 , 'Cutoff',  20, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  8 );
    ERP = pop_savemyerp(ERP, 'erpname', [SUB{i} '_audi_erp_ar_diff_waves_lpfilt'], 'filename', fullfile(Subject_Path, [SUB{i} '_audi_erp_ar_diff_waves_lpfilt.erp']));
    
%End subject loop
end 
    
%**********************************************************************************************************************************************************************
