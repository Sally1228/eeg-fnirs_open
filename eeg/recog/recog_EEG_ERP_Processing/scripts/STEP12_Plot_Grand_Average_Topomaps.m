%Script #12
%Operates on data averaged across participants
%Uses the output from Script #10: Grand_Average_ERPs.m
%This script loads the grand average ERP waveforms from Script #9, plots the topographic maps of the mean amplitude for the parent waveforms and difference waveform during the time window of the component, 
%and saves pdfs of all of the plots in the grand average ERPs folder.

%% Config
close all; clearvars;

% Assumptions & Notes:
% - Grand average ERP outputs are stored in outputs/Grand_Average_ERPs.
% - Channel location file is stored in inputs/standard-10-5-cap385.elp.

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);

input_dir = fullfile(project_dir, 'inputs');
output_dir = fullfile(project_dir, 'outputs');
ga_output_dir = fullfile(output_dir, 'Grand_Average_ERPs');


fprintf('[%s] Assumption: GA outputs in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ga_output_dir);
fprintf('[%s] Assumption: channel locations in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), fullfile(input_dir, 'standard-10-5-cap385.elp'));
fprintf('[%s] STEP12 Plot grand average topomaps\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%*************************************************************************************************************************************

%Set mean amplitude time window in milliseconds (e.g., 300 to 500 ms) 
timewindow = [300 500];

%Set EEG channels to include in the topomaps
chans = [1:60];

%Set bins to create corresponding topomaps; a separate topomap will be created for each bin specified
bins = [1 2 3]; 

%Set color scale limits for the topomaps for each bin in microvolts. To have scale limits set automatically, replace the code with myclim  = [];
myclim = [-2 4; -2 4; -2 4];

%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Load the grand average ERP waveforms outputted from Script #9 in ERPLAB format
ERP = pop_loaderp('filename', 'GA_recog_erp_ar_diff_waves.erp', 'filepath', ga_output_dir);  

%Look up channel location 
ERP = erpchanedit(ERP, fullfile(input_dir, 'standard-10-5-cap385.elp'));

%Convert mean amplitude time window to corresponding data points
[~,t_temp,~] = closest(ERP.times,timewindow);
timerange = t_temp(1):t_temp(2);

%Create figure
figure('Position',[500 500 1500 500]);

%Calculate the mean amplitude for each bin and channel and create a topomap of the mean amplitude for each bin
for b = 1:length(bins) 
    data2plot = squeeze(mean(ERP.bindata(chans,timerange,bins(b)),2)); 
    if isempty(myclim)
        clim(b,:) = round([min(data2plot) max(data2plot)]/0.5)*0.5; 
    else 
        clim = myclim; 
    end
    
    subplot(1,length(bins),b)
    topoplot(data2plot,ERP.chanlocs(chans),'maplimits',clim(b,:),'colormap',jet); 
    c = colorbar;  
    set(c,'YLim',clim(b,:),'fontsize',12);
    text(0, -0.8, ERP.bindescr{bins(b)},'fontsize',20, 'VerticalAlignment','bottom', ...
  'HorizontalAlignment', 'center'); 
end

%Save a pdf of the topomaps
save2pdf(fullfile(ga_output_dir, 'GA_recog_Topomaps.pdf'));
close all

%*************************************************************************************************************************************
