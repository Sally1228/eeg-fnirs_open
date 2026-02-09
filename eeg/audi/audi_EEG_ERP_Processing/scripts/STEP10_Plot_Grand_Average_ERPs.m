%Script #10
%Operates on data averaged across participants
%Uses the output from Script #9: Grand_Average_ERPs.m
%This script loads the low-pass filtered grand average ERP waveforms from Script #9, plots the difference waveforms, parent waveforms, ICA-corrected and uncorrected HEOG, and ICA-corrected VEOG, 
%and saves pdfs of all of the plots in the grand average ERPs folder.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Grand average ERP outputs are stored in outputs/Grand_Average_ERPs.

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'outputs');
ga_output_dir = fullfile(output_dir, 'Grand_Average_ERPs');


fprintf('[%s] Assumption: GA outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ga_output_dir);
fprintf('[%s] STEP10 Plot grand average ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


%List of subjects to process, based on the name of the folder that contains that subject's data
%sub-010质量不好
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-006', 'sub-008', 'sub-009'};   
task ='audi';

%*************************************************************************************************************************************

%Set baseline correction period in milliseconds
baselinecorr = '-200 0';

%Set x-axis scale in milliseconds
xscale = [-200.0 800.0   -200:200:800];

%Set y-axis scale in microvolts for the EEG channels for the parent waves
yscale_EEG_parent = [-10.0 15.0   -10:5:15];

%Set y-axis scale in microvolts for the EEG channels for the difference waves
yscale_EEG_diff = [-15.0 10.0   -15:5:10];

%Set y-axis scale in microvolts for the ICA-corrected and uncorrected bipolar HEOG channels
yscale_HEOG = [-15.0 15.0   -15:5:15];

%Set y-axis scale in microvolts for the ICA-corrected monopolar VEOG signals and corrected bipolar VEOG signal
yscale_VEOG = [-25.0 25.0   -25:10:25];

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Load the low-pass filtered grand average ERP waveforms outputted from Script #9 in .erp ERPLAB file format
ERP = pop_loaderp('filename', 'GA_audi_erp_ar_diff_waves_lpfilt.erp', 'filepath', ga_output_dir);    

%Plot the P3 rare and frequent parent waveforms at the key electrode sites of interest (FCz= 26 Cz=35, CPz=44, Pz=53)
ERP = pop_ploterps( ERP, [1 2], [26 35 44 53] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
save2pdf(fullfile(ga_output_dir, 'GA_audi_Parent_Waves.pdf'));
close all

%Plot the P3 rare-minus-frequent difference waveform at the key electrode sites of interest (FCz, Cz, CPz, Pz)
ERP = pop_ploterps( ERP, [3], [26 35 44 53] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_audi_Difference_Wave.pdf'));
close all

%Plot the P3 rare-minus-frequent difference waveform at the key electrode sites of interest (FCz, Cz, CPz, Pz) with the standard error of the mean (SEM)
ERP = pop_ploterps( ERP, [3], [26 35 44 53] , 'SEM', 'on', 'Transparency',  0.8, 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_audi_Difference_Wave_SEM.pdf'));
close all

%Plot the P3 rare and frequent parent waveforms at all electrode sites
ERP = pop_ploterps( ERP, [1 2], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
save2pdf(fullfile(ga_output_dir, 'GA_audi_Parent_Waves_All_Channels.pdf'));
close all

%Plot the P3 rare-minus-frequent difference waveform at all electrode sites
ERP = pop_ploterps( ERP, [3], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_audi_Difference_Wave_All_Channels.pdf'));
close all

% %Plot the parent (rare and frequent conditions) ICA-corrected and uncorrected HEOG signals 
% ERP = pop_ploterps( ERP, [1 2], [32 34] , 'Box', [1 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_HEOG);
% save2pdf([ga_output_dir filesep 'GA_P3_HEOG.pdf']);
% close all
% 
% %Plot the parent (rare and frequent conditions) ICA-corrected monopolar VEOG signals and corrected bipolar VEOG signal
% ERP = pop_ploterps( ERP, [1 2], [15 31 33] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_VEOG);
% save2pdf([ga_output_dir filesep 'GA_P3_VEOG.pdf']);
% close all

%*************************************************************************************************************************************
