%% STEP2_qc.m
% QC visualization for fNIRS preprocessing steps (script-only)
%
% What this script shows:
%   1) Raw intensity (by wavelength, mean across channels)
%   2) Raw OD (by wavelength, mean across channels)
%   3) Processed OD (after motion correction + filter)
%   4) HbO/HbR/HbT (mean across channels)
%   5) Motion masks (tIncAuto / tIncAutoCh) and channel mask (mlActAuto)
%   6) SCI summary (if available in qc)
%
% Outputs saved to: rootDir/derivatives/figures/

%% ===================== [1] CONFIG =====================
clear; clc;
rng(0);

% (A) Homer3 source path (edit this if needed)
% If you already addpath(genpath(Homer3)) elsewhere, set homer3Dir = ''.
homer3Dir = '';

% (B) Data root derived from this script's location
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
rootDir = fileparts(thisDir);

% (C) Output directory for QC figures (per-subject)
outDir = '';

% (D) Plot options
cfg = struct();
cfg.maxPlotWl = 3;
cfg.maxPlotCh = 20;    % for channel mask visualization
cfg.useLog10Intensity = true;

%% ===================== [2] ADD HOMER3 TO PATH (OPTIONAL) =====================
if ~isempty(homer3Dir)
    if exist(homer3Dir, 'dir') ~= 7
        error('homer3Dir does not exist: %s', homer3Dir);
    end
    addpath(genpath(homer3Dir));
end

%% ===================== [3] DISCOVER FILES =====================
subDirs = dir(fullfile(rootDir, 'sub-*'));
if isempty(subDirs)
    error('No sub-* folders found under rootDir: %s', rootDir);
end

%% ===================== [4] MAIN LOOP =====================
for si = 1:numel(subDirs)
    subName = subDirs(si).name;
    subPath = fullfile(rootDir, subName);
    outDir = fullfile(subPath, 'figures');
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    tmp = dir(fullfile(subPath, '*.snirf'));
    if isempty(tmp)
        fprintf('[SKIP] %s: no .snirf found\n', subName);
        continue;
    end

    for fi = 1:numel(tmp)
        inFile = fullfile(tmp(fi).folder, tmp(fi).name);
        [~, runBase, ~] = fileparts(inFile);
        runTag = sprintf('%s_%s', subName, runBase);

        outMat = fullfile(subPath, sprintf('%s_task-motor_fnirs_preproc.mat', subName));
        if exist(outMat, 'file') ~= 2
            fprintf('[SKIP] %s: preproc MAT not found: %s\n', runTag, outMat);
            continue;
        end

        fprintf('\n[QC] %s\n', runTag);

        %% ---- Load SNIRF raw intensity ----
        snirfObj = SnirfClass(inFile);
        if isempty(snirfObj.data)
            warning('SNIRF has empty data: %s', inFile);
            continue;
        end
        intensity = snirfObj.data(1);
        t = double(intensity.time);
        d = double(intensity.dataTimeSeries);
        nT = size(d,1);

        % measurement list info
        src = []; det = []; wl = [];
        try
            ml = intensity.measurementList;
            if ~isempty(ml)
                if isprop(ml(1),'sourceIndex'); src = [ml.sourceIndex]'; end
                if isprop(ml(1),'detectorIndex'); det = [ml.detectorIndex]'; end
                if isprop(ml(1),'wavelengthIndex'); wl = [ml.wavelengthIndex]'; end
            end
        catch
        end
        if isempty(wl) || numel(wl) ~= size(d,2)
            wl = ones(size(d,2),1);
        end
        wlVals = unique(wl);
        if numel(wlVals) > cfg.maxPlotWl
            wlVals = wlVals(1:cfg.maxPlotWl);
        end

        %% ---- Compute raw OD ----
        dod_raw = [];
        dod_raw_data = [];
        try
            dod_raw = hmrR_Intensity2OD(intensity);
            dod_raw_data = double(dod_raw.dataTimeSeries);
        catch
            warning('hmrR_Intensity2OD failed for %s', runTag);
        end

        %% ---- Load preproc MAT ----
        S = load(outMat);
        if isfield(S, 'dod')
            dod_proc = S.dod;
            dod_proc_data = double(dod_proc.dataTimeSeries);
        else
            dod_proc_data = [];
        end
        if isfield(S, 'dc')
            dc = S.dc;
            dc_data = double(dc.dataTimeSeries);
            dc_ml = [];
            try
                dc_ml = dc.measurementList;
            catch
            end
        else
            dc = [];
            dc_data = [];
            dc_ml = [];
        end
        if isfield(S, 'mlActAuto'); mlActAuto = S.mlActAuto; else; mlActAuto = []; end
        if isfield(S, 'tIncAuto'); tIncAuto = S.tIncAuto; else; tIncAuto = []; end
        if isfield(S, 'tIncAutoCh'); tIncAutoCh = S.tIncAutoCh; else; tIncAutoCh = []; end
        if isfield(S, 'qc'); qc = S.qc; else; qc = struct(); end

        %% ===================== FIG 1: OVERVIEW =====================
        fig1 = figure('Color','w','Position',[100 100 1100 700]);
        subplot(2,2,1); hold on;
        for wi = 1:numel(wlVals)
            idx = (wl == wlVals(wi));
            y = mean(d(:,idx), 2, 'omitnan');
            if cfg.useLog10Intensity
                y = log10(max(y, eps));
                ylabel('log10(Intensity)');
            else
                ylabel('Intensity');
            end
            plot(t, y, 'LineWidth', 1);
        end
        title('Raw intensity (mean by wavelength)');
        xlabel('Time (s)');
        grid on;
        legend(arrayfun(@(x)sprintf('wl%d',x), wlVals, 'UniformOutput', false), 'Location','best');

        subplot(2,2,2); hold on;
        if ~isempty(dod_raw_data)
            for wi = 1:numel(wlVals)
                idx = (wl == wlVals(wi));
                y = mean(dod_raw_data(:,idx), 2, 'omitnan');
                plot(t, y, 'LineWidth', 1);
            end
        end
        title('Raw OD (mean by wavelength)');
        xlabel('Time (s)');
        ylabel('OD');
        grid on;

        subplot(2,2,3); hold on;
        if ~isempty(dod_proc_data)
            for wi = 1:numel(wlVals)
                idx = (wl == wlVals(wi));
                if size(dod_proc_data,2) == numel(wl)
                    y = mean(dod_proc_data(:,idx), 2, 'omitnan');
                else
                    y = mean(dod_proc_data, 2, 'omitnan');
                end
                plot(t, y, 'LineWidth', 1);
            end
        end
        title('Processed OD (after TDDR+wavelet+filter)');
        xlabel('Time (s)');
        ylabel('OD');
        grid on;

        subplot(2,2,4); hold on;
        if ~isempty(dc_data)
            dataTypeLabel = [];
            try
                if ~isempty(dc_ml) && isprop(dc_ml(1), 'dataTypeLabel')
                    dataTypeLabel = {dc_ml.dataTypeLabel}';
                end
            catch
            end
            if ~isempty(dataTypeLabel) && numel(dataTypeLabel) == size(dc_data,2)
                prefLabels = {'HbO','HbR','HbT'};
                uniqLabels = {};
                for i = 1:numel(dataTypeLabel)
                    if ~isempty(dataTypeLabel{i}) && ~any(strcmp(uniqLabels, dataTypeLabel{i}))
                        uniqLabels{end+1} = dataTypeLabel{i}; %#ok<AGROW>
                    end
                end
                plotLabels = [prefLabels, uniqLabels];
                keep = false(size(plotLabels));
                for i = 1:numel(plotLabels)
                    if any(strcmp(dataTypeLabel, plotLabels{i}))
                        keep(i) = true;
                    end
                end
                plotLabels = plotLabels(keep);
                plotLabels = unique(plotLabels, 'stable');
                for i = 1:numel(plotLabels)
                    idx = strcmp(dataTypeLabel, plotLabels{i});
                    if any(idx)
                        y = mean(dc_data(:,idx), 2, 'omitnan');
                        plot(t, y, 'LineWidth', 1);
                    end
                end
                legend(plotLabels, 'Location','best');
            else
                y = mean(dc_data, 2, 'omitnan');
                plot(t, y, 'LineWidth', 1);
            end
        end
        title('Hb (mean across channels)');
        xlabel('Time (s)');
        ylabel('Conc (a.u.)');
        grid on;

        sgtitle(sprintf('%s | Overview', runTag), 'Interpreter','none');
        fig1Path = fullfile(outDir, sprintf('%s_step-qc_overview.png', runBase));
        saveas(fig1, fig1Path);
        close(fig1);

        %% ===================== FIG 2: QC MASKS =====================
        fig2 = figure('Color','w','Position',[100 100 1100 700]);
        subplot(2,2,1);
        mlActVec = [];
        if iscell(mlActAuto) && ~isempty(mlActAuto)
            mlActVec = mlActAuto{1};
        elseif ~isempty(mlActAuto)
            mlActVec = mlActAuto;
        end
        if ~isempty(mlActVec)
            plot(mlActVec, 'k.');
            ylim([-0.1 1.1]);
            ylabel('Active (1/0)');
            title('Channel activity mask (mlActAuto)');
        else
            text(0.1, 0.5, 'mlActAuto not found', 'FontSize', 12);
            axis off;
        end
        xlabel('Measurement index');
        grid on;

        subplot(2,2,2);
        tIncVec = [];
        if iscell(tIncAuto) && ~isempty(tIncAuto)
            tIncVec = tIncAuto{1};
        elseif ~isempty(tIncAuto)
            tIncVec = tIncAuto;
        end
        if ~isempty(tIncVec)
            plot(t, double(tIncVec), 'k-', 'LineWidth', 1);
            ylim([-0.1 1.1]);
            ylabel('tIncAuto');
            title('Motion mask (time)');
            grid on;
        else
            text(0.1, 0.5, 'tIncAuto not found', 'FontSize', 12);
            axis off;
        end
        xlabel('Time (s)');

        subplot(2,2,3);
        tIncCh = [];
        if iscell(tIncAutoCh) && ~isempty(tIncAutoCh)
            tIncCh = tIncAutoCh{1};
        elseif ~isempty(tIncAutoCh)
            tIncCh = tIncAutoCh;
        end
        if ~isempty(tIncCh)
            if size(tIncCh,1) ~= nT && size(tIncCh,2) == nT
                tIncCh = tIncCh';
            end
            tIncCh = double(tIncCh > 0);
            imagesc(tIncCh');
            colormap(gca, gray);
            caxis([0 1]);
            xlabel('Time (samples)');
            ylabel('Channel');
            title('Motion mask by channel (tIncAutoCh)');
        else
            text(0.1, 0.5, 'tIncAutoCh not found', 'FontSize', 12);
            axis off;
        end

        subplot(2,2,4);
        if isfield(qc, 'flag_motionDetectFailed')
            txt = sprintf('motionDetectFailed: %d\\ncritical: %d', ...
                qc.flag_motionDetectFailed, qc.flag_motionDetectCritical);
        else
            txt = 'motionDetect flags not found';
        end
        text(0.1, 0.6, txt, 'FontSize', 12);
        axis off;
        title('QC flags (motion)');

        sgtitle(sprintf('%s | QC masks', runTag), 'Interpreter','none');
        fig2Path = fullfile(outDir, sprintf('%s_step-qc_masks.png', runBase));
        saveas(fig2, fig2Path);
        close(fig2);

        %% ===================== FIG 3: SCI SUMMARY =====================
        if isfield(qc, 'sciMed') && ~isempty(qc.sciMed)
            fig3 = figure('Color','w','Position',[100 100 900 400]);
            sciMed = qc.sciMed(:);
            bar(sciMed);
            hold on;
            if isfield(qc, 'sciThr')
                yline(qc.sciThr, 'r--', 'LineWidth', 1.5);
            end
            xlabel('Src-Det pair index');
            ylabel('SCI (median)');
            title('SCI by src-det pair');
            grid on;
            fig3Path = fullfile(outDir, sprintf('%s_step-qc_sci.png', runBase));
            saveas(fig3, fig3Path);
            close(fig3);
        end
    end
end
