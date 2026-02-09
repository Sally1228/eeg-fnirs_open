%Script #9
%Operates on individual subject data
%Uses the output from Script #8: Average_ERPs.m
%This script loads the low-pass filtered averaged ERP waveforms from Script #7, plots the difference waveforms, parent waveforms, ICA-corrected and uncorrected HEOG, and ICA-corrected VEOG,
%and saves pdfs of all of the plots in the graphs folder located within each subjects's data folder.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Subject folders are stored in the sibling data directory (eeg/recog/sub-XXX).

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);
data_dir = fileparts(project_dir);


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP9 Plot individual ERPs\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = { 'sub-002', 'sub-003', 'sub-004', 'sub-005', 'sub-008'};    
task ='recog';

%**********************************************************************************************************************************************************************

%Set baseline correction period in milliseconds
baselinecorr = '-200 0';

%Set x-axis scale in milliseconds
xscale = [-200.0 1200.0   -200:200:1200];

%Set y-axis scale in microvolts for the EEG channels for the parent waves
yscale_EEG_parent = [-20.0 50.0   -20:10:50];

%Set y-axis scale in microvolts for the EEG channels for the difference waves
yscale_EEG_diff = [-20.0 30.0   -20:10:30];

%Set y-axis scale in microvolts for the ICA-corrected and uncorrected bipolar HEOG channels
yscale_HEOG = [-15.0 15.0   -15:5:15];

%Set y-axis scale in microvolts for the ICA-corrected monopolar VEOG signals and corrected bipolar VEOG signal
yscale_VEOG = [-25.0 25.0   -25:10:25];

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});
    graphs_dir = fullfile(Subject_Path, 'graphs');
    if ~exist(graphs_dir, 'dir')
        mkdir(graphs_dir);
    end

    %Load the low-pass filtered averaged ERP waveforms outputted from Script #7 in .erp ERPLAB file format
    ERP = pop_loaderp('filename', [SUB{i} '_recog_erp_ar_diff_waves_lpfilt.erp'], 'filepath', Subject_Path);    
    
    %Plot the recog FN400 OldHit and NewReject parent waveforms at the key electrode sites of interest (Fz=15 /FCz= 26/Cz=35)
    ERP = pop_ploterps( ERP, [1 2], [15 26 35] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
    save2pdf(fullfile(graphs_dir, [SUB{i} '_recog_FN400_Parent_Waves.pdf']));
    close all

    %Plot the PFN400 NewRejectt-OldHit difference waveform at the key electrode sites of interest (FCz= 26 Cz=35, CPz=44, Pz=53)
    %[3] 表示只画 第 3 个 bin 的波形；[15 26 35] 只画这几个通道；'Box', [2 2] 2×2排版；'blc', baselinecorr 基线校正
    ERP = pop_ploterps( ERP, [3], [15 26 35] , 'Box', [2 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
    save2pdf(fullfile(graphs_dir, [SUB{i} '_recog_FN400_Difference_Wave.pdf']));
    close all

%%     可以选是否画所有channel 
%     %Plot the FN400 OldHit and NewReject parent waveforms at all electrode sites
%     ERP = pop_ploterps( ERP, [1 2], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
%     save2pdf([Subject_Path 'graphs' filesep SUB{i} '_recog_FN400_Parent_Waves_All_Channels.pdf']);
%     close all
% 
%     %Plot the PFN400 NewRejectt-OldHit difference waveform at all electrode sites
%     ERP = pop_ploterps( ERP, [3], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
%     save2pdf([Subject_Path 'graphs' filesep SUB{i} '_recog_FN400_Difference_Wave_All_Channels.pdf']);
%     close all

    %Plot the recog LPC OldHit and NewReject parent waveforms at the key electrode sites of interest (P3=47 /Pz=53 /P4 = 48 + CP3=38/CPz=44/CP4=39)
    ERP = pop_ploterps( ERP, [1 2], [ 47 53 48 38 44 39 ] , 'Box', [3 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
    save2pdf(fullfile(graphs_dir, [SUB{i} '_recog_LPC_Parent_Waves.pdf']));
    close all

    %Plot the LPC NewRejectt-OldHit difference waveform at the key electrode sites of interest (FCz= 26 Cz=35, CPz=44, Pz=53)
    %[3] 表示只画 第 3 个 bin 的波形；[ 47 53 48 38 44 39] 只画这几个通道；'Box', [2 3] 2×3排版；'blc', baselinecorr 基线校正
    ERP = pop_ploterps( ERP, [3],  [ 47 53 48 38 44 39 ] , 'Box', [3 2], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
    save2pdf(fullfile(graphs_dir, [SUB{i} '_recog_LPC_Difference_Wave.pdf']));
    close all

%%     可以选是否画所有channel 
%     %Plot the LPC OldHit and NewReject parent waveforms at all electrode sites
%     ERP = pop_ploterps( ERP, [1 2], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_parent);
%     save2pdf([Subject_Path 'graphs' filesep SUB{i} '_recog_LPC_Parent_Waves_All_Channels.pdf' ]);
%     close all
% 
%     %Plot the LPC NewRejectt-OldHit difference waveform at all electrode sites
%     ERP = pop_ploterps( ERP, [3], [1:60] , 'Box', [8 8], 'blc', baselinecorr, 'Maximize', 'on', 'Style', 'Classic', 'xscale', xscale,  'yscale', yscale_EEG_diff);
%     save2pdf([Subject_Path 'graphs' filesep SUB{i} '_recog_LPC_Difference_Wave_All_Channels.pdf']);
%     close all


%End subject loop
end

%*************************************************************************************************************************************
