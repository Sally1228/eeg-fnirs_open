%Script #11
%Operates on data averaged across participants
%Uses the output from Script #10: Grand_Average_ERPs.m
%This script loads the low-pass filtered grand average ERP waveforms from Script #9, plots the difference waveforms, parent waveforms, ICA-corrected and uncorrected HEOG, and ICA-corrected VEOG, 
%and saves pdfs of all of the plots in the grand average ERPs folder.


%% Config
close all; clearvars;

% Assumptions & Notes:
% - Grand average ERP outputs are stored in outputs/Grand_Average_ERPs.

script_dir = fileparts(mfilename('fullpath'));

project_dir= fileparts(script_dir);

output_dir = fullfile(project_dir, 'outputs');
ga_output_dir = fullfile(output_dir, 'Grand_Average_ERPs');


fprintf('[%s] Assumption: GA outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ga_output_dir);
fprintf('[%s] STEP11 Plot grand average ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-008'};    
task ='recog';

%*************************************************************************************************************************************

%Set baseline correction period in milliseconds
baselinecorr = '-200 0';

%Set x-axis scale in milliseconds
xscale = [-200.0 1200.0   -200:200:1200];

%Set y-axis scale in microvolts for the EEG channels for the parent waves
yscale_EEG_parent = [-10.0 15.0   -10:5:15];

%Set y-axis scale in microvolts for the EEG channels for the difference waves
yscale_EEG_diff = [-15.0 10.0   -15:5:10];


%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Load the low-pass filtered grand average ERP waveforms outputted from Script #9 in .erp ERPLAB file format
ERP = pop_loaderp('filename', 'GA_recog_erp_ar_diff_waves_lpfilt.erp', 'filepath', ga_output_dir);    

%Plot the FN400 parent waveforms at the key electrode sites of interest (Fz=15 /FCz= 26/Cz=35)
ERP = pop_ploterps( ERP, [1 2], [15 26 35] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
save2pdf(fullfile(ga_output_dir, 'GA_recog_FN400_Parent_Waves.pdf'));
close all
%Plot the FN400 difference waveform at the key electrode sites of interest (Fz=15 /FCz= 26/Cz=35)
ERP = pop_ploterps( ERP, [3], [15 26 35] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_recog_FN400_Difference_Wave.pdf'));
close all
%Plot the FN400 difference waveform at the key electrode sites of interest (FCz, Cz, CPz, Pz) with the standard error of the mean (SEM)
ERP = pop_ploterps( ERP, [3], [15 26 35] , 'SEM', 'on', 'Transparency',  0.8, 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_recog_FN400_Difference_Wave_SEM.pdf'));
close all

% % 可以选择是否画所有导联 Plot the FN400 parent waveforms at all electrode sites
% ERP = pop_ploterps( ERP, [1 2], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
% save2pdf([ga_output_dir filesep 'GA_recog_FN400_Parent_Waves_All_Channels.pdf']);
% close all
% %Plot the FN400 difference waveform at all electrode sites
% ERP = pop_ploterps( ERP, [3], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
% save2pdf([ga_output_dir filesep 'GA_recog_FN400_Difference_Wave_All_Channels.pdf']);
% close all

%Plot the LPC parent waveforms at the key electrode sites of interest (P3=47 /Pz=53 /P4 = 48 + CP3=38/CPz=44/CP4=39)
ERP = pop_ploterps( ERP, [1 2], [ 47 53 48 38 44 39 ] , 'Box', [3 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
save2pdf(fullfile(ga_output_dir, 'GA_recog_LPC_Parent_Waves.pdf'));
close all
%Plot the LPC difference waveform at the key electrode sites of interest (P3=47 /Pz=53 /P4 = 48 + CP3=38/CPz=44/CP4=39)
ERP = pop_ploterps( ERP, [3], [ 47 53 48 38 44 39 ] , 'Box', [3 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_recog_LPC_Difference_Wave.pdf'));
close all
%Plot the LPC difference waveform at the key electrode sites of interest (FCz, Cz, CPz, Pz) with the standard error of the mean (SEM)
ERP = pop_ploterps( ERP, [3], [ 47 53 48 38 44 39 ] , 'SEM', 'on', 'Transparency',  0.8, 'Box', [3 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
save2pdf(fullfile(ga_output_dir, 'GA_recog_LPC_Difference_Wave_SEM.pdf'));
close all

% % 可以选择是否画所有导联 Plot the LPC parent waveforms at all electrode sites
% ERP = pop_ploterps( ERP, [1 2], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
% save2pdf([ga_output_dir filesep 'GA_recog_LPC_Parent_Waves_All_Channels.pdf']);
% close all
% %Plot the LPC difference waveform at all electrode sites
% ERP = pop_ploterps( ERP, [3], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
% save2pdf([ga_output_dir filesep 'GA_recog_LPC_Difference_Wave_All_Channels.pdf']);
% close all

close all

%*************************************************************************************************************************************
