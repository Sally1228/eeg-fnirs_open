%Script #13 
%Operates on individual subject data
%Uses the output from Script #12: Measure_ERPs.m
%This script loads the individual subject difference waveform measurement values from Script #12, creates histogram plots of the single-participant measurement values for each measure of amplitude and latency (e.g., mean amplitude, peak amplitude), 
%and saves pdfs of all of the plots in the ERP Measurements folder.


%% Config
close all; clearvars;

script_dir = fileparts(mfilename('fullpath'));

project_dir = fileparts(script_dir);
output_dir = fullfile(project_dir, 'outputs');
erp_output_dir = fullfile(output_dir, 'ERP_Measurements');


fprintf('[%s] Assumption: ERP measurements in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), erp_output_dir);
fprintf('[%s] STEP13 Plot measurement histograms\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%*************************************************************************************************************************************

%Specify the measurements that will be used to create histograms using the text files created by Script #12 (e.g., Mean_Amplitude_Diff_Waves.txt)
Measures = {'Mean_Amplitude','Peak_Amplitude','Peak_Latency','50%_Area_Latency','Onset_Latency'}; 

%Set the histogram bin edges (lower and upper limits for each bar) to be plotted for each measurement
Xranges = {0:0.5:4, 3:1:9, 300:50:500, 300:20:400, 150:50:350};

%Loop through each measurement listed in Measures
for m = 1:length(Measures)
    
    meas = Measures{m};
    
    %Calculate center points of histogram bins based on bin edges specified above
    intervals = diff(Xranges{m});
    bincenters = Xranges{m}(1:end-1) + intervals/2;
    
    %Load the single-participant measurement values outputted from Script #12
    Meas_file = fullfile(erp_output_dir, [meas '_Diff_Waves_audi.txt']);
    Meas_table = readtable(Meas_file, 'delimiter', '\t', 'HeaderLines', 1);
    Meas_data = Meas_table{:,1};
    
    %Create histogram plots
    [bincounts] = histcounts(Meas_data,Xranges{m});
    bar(bincenters,bincounts,0.75,'FaceColor',[0.2 0.2 0.5]);
    ylim([0 30]); xticks(bincenters); 
    ylabel('Number of Subjects')
    xlabel(strrep(meas,'_',' '))
    set(gca,'fontsize',12)
    legend({'Rare-Frequent'});
    
    %Saves the histogram plots to pdf
    save2pdf(fullfile(erp_output_dir, [meas '_Histogram_audi.pdf']));
    close all
    
%End measures loop    
end

%*************************************************************************************************************************************
