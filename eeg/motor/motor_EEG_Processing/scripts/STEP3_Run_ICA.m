% %Script #3
% %Operates on individual subject data
% %Uses the output from Script #2: ICA_Prep.m
% %This script loads the outputted semi-continuous EEG data file from Script #2, computes the ICA weights that will be used for artifact correction of ocular artifacts, transfers the ICA weights to the 
% %continuous EEG data file outputted from Script #1 (e.g., without the break periods and noisy segments of EEG removed), and saves a pdf of the topographic maps of the ICA weights.
%

%% Config
close all; clearvars;

% script_dir =  '/Users/zhaifeifei/Desktop/eeg_fnirs/eeg/motor/motor_EEG_ERP_Processing/scripts'
script_dir = fileparts(mfilename('fullpath'));

project_dir= fileparts(script_dir);
input_dir = fullfile(project_dir, 'inputs');
data_dir = fileparts(project_dir);


fprintf('[%s] Assumption: subject folders live in %s\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), data_dir);
fprintf('[%s] STEP3 Run ICA\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'sub-001','sub-002','sub-003', 'sub-004', 'sub-005', 'sub-006','sub-007', 'sub-008', 'sub-009', 'sub-010' };
task ='motor';

%***********************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    fprintf('[%s] Processing %s (%d/%d)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), SUB{i}, i, length(SUB));
    
    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = fullfile(data_dir, SUB{i});

    %Load the semi-continuous EEG data file outputted from Script #2 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_pre2.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_pre2'], 'overwrite','on', 'gui', 'off');   

    %Compute ICA weights with binICA (a compiled and faster version of ICA). If binICA is not an option (e.g., on a Windows machine), use runICA by replacing the code with the following: 
    %对1-60通道跑ica，没有包括耳电极
    EEG = pop_runica(EEG,'extended',1,'chanind', [1:60]);
    %Note that the bipolar HEOG and VEOG channels are not included in the channel list for computing ICA weights, because they are not linearly independent of the channels that were used to create them
    %用下面这个会报错，没有binica，用上面
    % EEG = pop_runica(EEG,'extended',1,'icatype','binica','chanind', [1:31]); 
   
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_pre2_weighted'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_pre2_weighted.set']), 'overwrite', 'on', 'gui', 'off');

    %Load the continuous EEG data file outputted from Script #1 in .set EEGLAB file format
    EEG = pop_loadset('filename', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_notch.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i}  '_task-' task  '_eeg_ds_reref_hpfilt_notch'], 'gui', 'off'); 

    %Transfer ICA weights to the continuous EEG data file (e.g., without the break periods and noisy segments of data removed)
    EEG = pop_editset(EEG, 'icachansind', 'ALLEEG(2).icachansind', 'icaweights', 'ALLEEG(2).icaweights', 'icasphere', 'ALLEEG(2).icasphere');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted'], ...
        'savenew', fullfile(Subject_Path, [SUB{i} '_task-' task '_eeg_ds_reref_hpfilt_ica_weighted.set']), 'overwrite', 'on', 'gui', 'off');

    %Run ICLabel if available (use current EEG with ICA weights)
    has_iclabel = exist('pop_iclabel','file') == 2 || exist('iclabel','file') == 2;
    if has_iclabel
        try
            EEG.srate = round(EEG.srate);
            EEG = eeg_checkset(EEG);
            EEG = pop_iclabel(EEG, 'default');
        catch ME
            warning('ICLabel failed: %s', ME.message);
            has_iclabel = false;
        end
        if ~(isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification'))
            has_iclabel = false;
        end
    else
        warning('ICLabel plugin not found; skipping ICLabel classification.');
    end

    %Save a pdf of the topographic maps of the ICA weights for later review
    set(groot,'DefaultFigureColormap',jet)
    graphs_dir = fullfile(Subject_Path, 'graphs');
    if ~exist(graphs_dir, 'dir')
        mkdir(graphs_dir);
    end

    comps = 1:size(EEG.icawinv, 2);
    nrows = ceil(sqrt(numel(comps)));
    ncols = ceil(numel(comps) / nrows);

    has_labels = has_iclabel && isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') ...
            && isfield(EEG.etc.ic_classification, 'ICLabel');
    if has_labels
        classes = EEG.etc.ic_classification.ICLabel.classes;
        if isstring(classes)
            classes = cellstr(classes);
        end
        class_probs = EEG.etc.ic_classification.ICLabel.classifications;
        [max_prob, max_idx] = max(class_probs, [], 2);
    end

    fig = figure('Color','w', 'Visible','off');
    tile_size = 240;
    fig_w = max(1600, ncols * tile_size);
    fig_h = max(1200, nrows * tile_size);
    set(fig, 'Units', 'pixels', 'Position', [100 100 fig_w fig_h]);

    for ci = 1:numel(comps)
        c = comps(ci);
        subplot(nrows, ncols, ci);
        topoplot(EEG.icawinv(:, c), EEG.chanlocs(EEG.icachansind), 'electrodes','on');
        if has_labels && c <= numel(max_idx)
            title(sprintf('IC %d: %s %.2f', c, classes{max_idx(c)}, max_prob(c)), 'FontSize', 7);
        else
            title(sprintf('IC %d', c), 'FontSize', 7);
        end
    end

    out_file = fullfile(graphs_dir, [SUB{i} '_motor_ICA_Weights.png']);
    exportgraphics(fig, out_file, 'Resolution', 150);
    close(fig);

    
%End subject loop
end
% 
% %***********************************************************************************************************************************************
