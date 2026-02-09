%% fnirs_preproc_homer3_snirf_3WL_NOFUNC.m
% 这个版本加入 SCI，运动伪迹用 TDDR + wavelet
%  Robust fNIRS preprocessing + QC pipeline using Homer3 (script-only, no functions)
%
% Design goals:
%   (1) MUST run with SNIRF that has 3 wavelengths
%   (2) Be readable: very explicit comments + step-by-step logging
%   (3) Robust to Homer3 "cell vs non-cell" input conventions
%   (4) Produce derivatives: preproc MAT + QC figures + QC_summary.csv
%
% What this script does (high-level):
%   A. Find files: rootDir/sub-*/\*.snirf (non-recursive within each subject folder)
%   B. For each file:
%      1) Load SNIRF -> intensity (raw) + probe
%      2) Detect number of wavelengths nW (IMPORTANT for PPF length)
%      3) Channel QC (PruneChannels) -> mlActAuto (keep/reject each measurement)
%      4) Convert intensity -> OD
%      5) Detect motion segments on OD -> tIncAuto / tIncAutoCh (mask)
%      6) Motion correction: TDDR -> Wavelet (abort if wavelet fails)
%      7) Bandpass filter on OD (0.01–0.5 Hz by default)
%      8) Convert OD -> HbO/HbR/HbT via MBLL (OD2Conc), PPF length = nW
%      9) Save outputs + QC figures + QC table row
%
% IMPORTANT:
%   - You must have Homer3 source on MATLAB path (SnirfClass, hmrR_* functions).
%   - If wavelet motion correction errors, the script will abort.

clear; clc;

%% ===================== [1] USER SETTINGS =====================
% (A) Homer3 source path (edit this if needed)
% If you already addpath(genpath(Homer3)) elsewhere, set homer3Dir = ''.
homer3Dir = '';   

% (A2) TDDR source path (folder containing TDDR.m)
% Example: fullfile(getenv('HOME'),'Downloads','TDDR-master')
tddrDir = '';

% (B) Data root derived from this script's location
thisDir = fileparts(mfilename('fullpath'));

rootDir = fileparts(thisDir);

%TDDR source path (folder containing TDDR.m)
tddrDir = fullfile(thisDir,'TDDR');

% (C) Output directory for QC summary (use script folder)
outDir = thisDir;

% (D) Processing parameters (widely used defaults / common starting points)
params = struct();

% ---- PruneChannels (measurement QC) ----
% We use auto-dRange by default because different devices scale intensity differently.
params.useAutoDRange = true;
params.fixedDRange   = [1e4, 1e7];     % 设置默认参数，或者是用自己设置的参数 use ONLY if you know your intensity scaling matches
params.SNRthresh     = 2;              % Homer3 common default
params.SDrange       = [0.0, 45.0];    % mm, typical adult cap range

params.useSCI        = true;
params.SCIband       = [0.5, 2.5];     % Hz, typical cardiac band
params.SCIwinSec     = 10;             % seconds (<=0 means whole-recording)
params.SCIthr        = 0.6;            % threshold on median SCI
params.SCIwlPair     = [1, 2];         % wavelength indices for SCI

% ---- Motion detection on OD ----
% If motion detection marks too much: increase STDEVthresh or AMPthresh
% If motion detection marks too little: decrease STDEVthresh or AMPthresh
params.tMotion     = 0.5;              % sec
params.tMask       = 1.0;              % sec
params.STDEVthresh = 20.0;             % common practical start (tighter than 50)
params.AMPthresh   = 0.5;              % OD amplitude threshold (scale-dependent)

% ---- Motion correction ----
params.doTDDR      = true;
params.doWavelet   = true;
params.iqrWavelet  = 1.5;              % typical

% ---- Bandpass filter (OD stage) ----
params.hpf = 0.01;                     % 滤波范围Hz
params.lpf = 0.20;                     % Hz (try 0.20 if you want stronger physiological suppression)

% ---- OD -> Hb concentrations ----
% PPF MUST match #wavelengths. We store as scalar; script will replicate to length nW.
params.ppfScalar = 6.0;                % start with 6.0, so output is "molar * ppf"

% (E) QC flags for reporting (non-positive intensity may trigger skip)
qcRule = struct();
qcRule.minActiveMeasFrac = 0.60;       % <60% measurements survive -> flag 活跃通道比例下限。比如 0.60 表示“保留的测量通道 <60%”就记为 flag_lowActive=1。
qcRule.maxMotionFrac     = 0.20;       % >20% time is motion -> flag   运动伪迹比例上限。比如 0.20 表示“被判为运动的时间 >20%”就记为 flag_highMotion=1。
qcRule.maxNonPositiveMeasFrac = 0.05;  % >5% measurements have <=0 intensity -> skip run   非正强度通道比例上限。比如 0.05 表示“含≤0强度的测量通道 >5%”就直接跳过该 run，并标记 flag_nonPositiveSevere=1。

% (F) Debug options
debug = struct();
debug.processOnlyFirstFile = false;    % set true to test quickly

%% ===================== [2] ADD HOMER3 TO PATH (OPTIONAL) =====================
if ~isempty(homer3Dir)
    if exist(homer3Dir, 'dir') ~= 7
        error('homer3Dir does not exist: %s', homer3Dir);
    end
    addpath(genpath(homer3Dir));
end
if ~isempty(tddrDir)
    if exist(tddrDir, 'dir') ~= 7
        error('tddrDir does not exist: %s', tddrDir);
    end
    addpath(genpath(tddrDir));
end

%% ===================== [2b] CHECK WAVELET AVAILABILITY =====================
waveletAvailable = (exist('hmrR_MotionCorrectWavelet', 'file') == 2);
if params.doWavelet && ~waveletAvailable
    error('WaveletMissing:MotionCorrectWavelet', ...
        'Wavelet requested but hmrR_MotionCorrectWavelet is not available.');
end

%% ===================== [2c] CHECK TDDR AVAILABILITY =====================
tddrAvailable = (exist('TDDR', 'file') == 2);
if params.doTDDR && ~tddrAvailable
    error('TDDRMissing:ExternalFunction', ...
        'TDDR requested but TDDR.m is not available on the MATLAB path.');
end

%% ===================== [3] FIND SUBJECTS =====================
subDirs = dir(fullfile(rootDir, 'sub-*'));
if isempty(subDirs)
    error('No sub-* folders found under rootDir: %s', rootDir);
end

% We'll accumulate QC rows into a cell array, then write to CSV at the end.
qcHeader = {'sub','run','inFile','outMat','nW','durationSec','nMeasTotal','nMeasActive', ...
            'activeMeasFrac','motionFrac','nonPositiveMeasFrac','nonPositiveMeasCount', ...
            'nFailSNRMeas','nFailSDMeas','nFailSCIMeas','nFailSNRPair','nFailSDPair','nFailSCIPair', ...
            'flag_nonPositiveCritical','flag_nonPositiveSevere','flag_motionDetectFailed', ...
            'flag_motionDetectCritical','didTDDR','didWavelet','sciMedByPair','sciGoodPctByPair', ...
            'sciThr','sciWlPair','dRangeByWl','ppfUsed','addedOffset', ...
            'flag_lowActive','flag_highMotion','note'};
qcData = cell(0, numel(qcHeader));

%% ===================== [4] MAIN LOOP =====================
for si = 1:numel(subDirs)

    subName = subDirs(si).name;
    subPath = fullfile(rootDir, subName);
    figDir = fullfile(subPath, 'figures');
    if ~exist(figDir, 'dir'); mkdir(figDir); end

    fprintf('\n==================== SUBJECT: %s ====================\n', subName);

    % --- Search .snirf file ---
    tmp = dir(fullfile(subPath, '*.snirf'));
    fileList = cell(numel(tmp), 1);
    for k = 1:numel(tmp)
        fileList{k} = fullfile(tmp(k).folder, tmp(k).name);
    end

    if isempty(fileList)
        fprintf('[SKIP] %s: no .snirf found (you said you have SNIRF; please check paths)\n', subName);
        continue;
    end

    for fi = 1:numel(fileList)

        inFile = fileList{fi};
        [~, runBase, ~] = fileparts(inFile);
        runTag = sprintf('%s_%s', subName, runBase);

        fprintf('\n--- RUN: %s ---\nFile: %s\n', runTag, inFile);

        note = "";
        try
            %% ========== STEP 0: LOAD SNIRF ==========
            % SnirfClass loads .snirf and provides:
            %   snirfObj.data(1): DataClass with time + intensity time series
            %   snirfObj.probe  : ProbeClass with optode geometry + wavelengths
            snirfObj = SnirfClass(inFile);

            if isempty(snirfObj.data)
                error('SNIRF has empty "data". Cannot proceed.');
            end

            intensity = snirfObj.data(1);
            probe     = snirfObj.probe;

            t = double(intensity.time);
            d = double(intensity.dataTimeSeries);

            if size(d,1) ~= numel(t)
                error('Time length (%d) does not match data rows (%d).', numel(t), size(d,1));
            end

            nT = size(d,1);
            nMeas = size(d,2);

            % Basic duration
            durationSec = t(end) - t(1);

            %% ========== STEP 0b: DETECT #WAVELENGTHS + MEAS LIST ==========
            % We try multiple ways; different Homer3 versions store this differently.
            nW = [];
            src = [];
            det = [];
            wl  = [];

            % infer from measurementList (preferred)
            try
                ml = intensity.measurementList;
                if ~isempty(ml)
                    if isprop(ml(1),'sourceIndex')
                        src = [ml.sourceIndex]';
                    end
                    if isprop(ml(1),'detectorIndex')
                        det = [ml.detectorIndex]';
                    end
                    if isprop(ml(1),'wavelengthIndex')
                        wl = [ml.wavelengthIndex]';
                        nW = numel(unique(wl));
                    end
                end
            catch
                % ignore, fallback below
            end

            % fallback: probe wavelengths
            if isempty(nW)
                try
                    if isprop(probe,'wavelengths') && ~isempty(probe.wavelengths)
                        nW = numel(probe.wavelengths);
                    end
                catch
                end
            end

            if isempty(nW)
                nW = 3;
                note = note + "[nW fallback=3] ";
            end

            fprintf('Detected wavelengths (nW) = %d\n', nW);

            %% ========== STEP 1: NON-POSITIVE INTENSITY QC (CRITICAL) ==========
            % OD uses -log(I/I0). Non-positive values are invalid and must be flagged.
            addedOffset = 0;
            nonPositiveMeas = false(nMeas,1);
            nonPositiveMeasCount = 0;
            nonPositiveMeasFrac = 0;
            flag_nonPositiveCritical = false;
            flag_nonPositiveSevere = false;

            nonPositiveMask = (d <= 0);
            if any(nonPositiveMask(:))
                nonPositiveMeas = any(nonPositiveMask, 1)';
                nonPositiveMeasCount = sum(nonPositiveMeas);
                nonPositiveMeasFrac = nonPositiveMeasCount / max(nMeas,1);
                flag_nonPositiveCritical = true;
                note = note + "[nonPositiveIntensity] ";
                fprintf('Non-positive intensity detected: %d/%d measurements (%.1f%%)\n', ...
                    nonPositiveMeasCount, nMeas, 100*nonPositiveMeasFrac);

                % Hard-remove affected measurements to avoid invalid OD
                d(:, nonPositiveMeas) = NaN;
                intensity.dataTimeSeries = d;

                if nonPositiveMeasFrac > qcRule.maxNonPositiveMeasFrac
                    flag_nonPositiveSevere = true;
                    note = note + "[nonPositiveSevere->skip] ";
                    fprintf('[SKIP] %s: non-positive intensity exceeds threshold (%.1f%% > %.1f%%)\n', ...
                        runTag, 100*nonPositiveMeasFrac, 100*qcRule.maxNonPositiveMeasFrac);

                    qcRow = {subName, runBase, inFile, '', nW, durationSec, nMeas, NaN, ...
                             NaN, NaN, nonPositiveMeasFrac, nonPositiveMeasCount, ...
                             NaN, NaN, NaN, NaN, NaN, NaN, ...
                             flag_nonPositiveCritical, flag_nonPositiveSevere, NaN, NaN, ...
                             false, false, '', '', NaN, '', '', '', addedOffset, NaN, NaN, char(note)};
                    qcData = [qcData; qcRow]; %#ok<AGROW>
                    continue;
                end
            end

            % Optional: compute source-detector distance (use 2D positions only)
            % 对每个通道计算欧氏距离 norm(srcPos - detPos)，放入 dist(ii)
            dist = nan(nMeas,1);
            if ~isempty(src) && ~isempty(det) && numel(src) == nMeas
                srcPos = [];
                detPos = [];
                if isprop(probe,'sourcePos2D') && ~isempty(probe.sourcePos2D)
                    srcPos = probe.sourcePos2D;
                end
                if isprop(probe,'detectorPos2D') && ~isempty(probe.detectorPos2D)
                    detPos = probe.detectorPos2D;
                end
                if ~isempty(srcPos) && ~isempty(detPos)
                    for ii = 1:nMeas
                        sIdx = src(ii); dIdx = det(ii);
                        if sIdx<=size(srcPos,1) && dIdx<=size(detPos,1)
                            dist(ii) = norm(srcPos(sIdx,:) - detPos(dIdx,:));
                        end
                    end
                end
            end

            %% ========== STEP 2: PREPARE MANUAL MASKS (NONE BY DEFAULT) ==========
            % Homer3 functions often expect tIncMan/mlActMan as cell arrays (for multiple runs).
            % Keep vector versions for fallbacks where supported.
            tIncMan_cell  = {ones(nT,1)}; % Time Inclusion Manual，1=保留，0=剔除
            mlActMan_cell = {ones(nMeas,1)}; %Measurement List Active Manual，每个测量通道1=保留，0=剔除
            tIncMan_vec   = ones(nT,1);
            mlActMan_vec  = ones(nMeas,1);
            if flag_nonPositiveCritical
                mlActMan_vec(nonPositiveMeas) = 0;
                mlActMan_cell = {mlActMan_vec};
                note = note + "[nonPositiveMeasInactive] ";
            end

            %% ========== STEP 2b: SCI (SCALP COUPLING INDEX) ==========
            % 计算 sci
            sciOK = true(nMeas,1);
            sciPairs = [];
            sciMed = [];
            sciGoodPct = [];
            sciComputed = false;
            if params.useSCI
                if isempty(src) || isempty(det) || isempty(wl) || numel(wl) ~= nMeas
                    note = note + "[SCI skip no measurementList] ";
                    %需要 Signal Processing Toolbox 的 butter/filtfilt。
                elseif ~exist('butter', 'file') || ~exist('filtfilt', 'file')
                    error('SCIMissing:SignalToolbox', ...
                        'SCI requires butter/filtfilt (Signal Processing Toolbox).');
                else
                    dt = median(diff(t));
                    if ~isfinite(dt) || dt <= 0
                        error('SCIInvalidTime:NonPositiveDT', 'Invalid time vector for SCI.');
                    end
                    fs = 1 / dt;
                    nyq = fs / 2;
                    hi = min(params.SCIband(2), nyq * 0.99); %心搏带带通滤波params.SCIband 默认 [0.5, 2.5] Hz
                    lo = min(params.SCIband(1), hi * 0.8);
                    lo = max(lo, eps);
                    if hi <= 0.2 || lo >= hi
                        error('SCIBandInvalid:LowFS', 'Sampling rate too low for SCI bandpass. fs=%.3f Hz', fs);
                    end
                    [b,a] = butter(3, [lo hi] / nyq, 'bandpass');%带通滤波

                    if isempty(params.SCIwinSec) || params.SCIwinSec <= 0
                        edges = [t(1) t(end) + eps];
                    else
                        edges = t(1):params.SCIwinSec:t(end);
                        if edges(end) < t(end)
                            edges(end+1) = t(end) + eps; %#ok<AGROW>
                        end
                    end
                    W = numel(edges) - 1;
        
                    % 以src‑det 为单位计算 SCI。用两个波长计算sci
                    sciPairs = unique([src det], 'rows');
                    Np = size(sciPairs, 1);
                    sci_win = nan(Np, W);

                    minFiltLen = 3 * max(numel(a), numel(b));
                    if nT <= minFiltLen
                        note = note + "[SCI skip short data] ";
                    else
                        for pi = 1:Np
                            s = sciPairs(pi,1); dIdx = sciPairs(pi,2);
                            idx1 = find(src == s & det == dIdx & wl == params.SCIwlPair(1), 1);
                            idx2 = find(src == s & det == dIdx & wl == params.SCIwlPair(2), 1);
                            if isempty(idx1) || isempty(idx2)
                                continue;
                            end
                            x0 = d(:, idx1);
                            y0 = d(:, idx2);

                            x0(x0 <= 0) = NaN;
                            y0(y0 <= 0) = NaN;
                            if any(~isfinite(x0)) || any(~isfinite(y0))
                                continue;
                            end

                            sx = std(x0, [], 'omitnan');
                            sy = std(y0, [], 'omitnan');
                            if ~isfinite(sx) || ~isfinite(sy) || sx == 0 || sy == 0
                                continue;
                            end
                            x0 = (x0 - mean(x0, 'omitnan')) / sx;
                            y0 = (y0 - mean(y0, 'omitnan')) / sy;

                            xf = filtfilt(b, a, x0);
                            yf = filtfilt(b, a, y0);

                            for w = 1:W
                                m = (t >= edges(w)) & (t < edges(w+1));
                                if nnz(m) < 10
                                    continue;
                                end
                                r = corr(xf(m), yf(m), 'rows', 'complete');
                                if isfinite(r)
                                    sci_win(pi, w) = max(0, r);
                                end
                            end
                        end
                    end

                    sciMed = median(sci_win, 2, 'omitnan');
                    sciGoodPct = mean(sci_win >= params.SCIthr, 2, 'omitnan');
                    sciComputed = true;

                    sciOK_byPair = true(Np,1);
                    validSci = ~isnan(sciMed);
                    sciOK_byPair(validSci) = (sciMed(validSci) >= params.SCIthr);
                    for pi = 1:Np
                        idx = (src == sciPairs(pi,1)) & (det == sciPairs(pi,2));
                        sciOK(idx) = sciOK_byPair(pi);
                    end
                end
            end

            %% ========== STEP 3: AUTO dRange FOR PRUNE ==========
            % dRange is device dependent. Auto-estimation is safer than hard-coded [1e4,1e7].
            dRangeByWl = [];
            wlVals = [];
            if params.useAutoDRange
                dMean = mean(d, 1, 'omitnan')';
                % wavelength-specific ranges
                if ~isempty(wl) && numel(wl) == nMeas
                    wlVals = unique(wl);
                    dRangeByWl = nan(numel(wlVals), 2);
                    for wi = 1:numel(wlVals)
                        dMeanWl = dMean(wl == wlVals(wi));
                        dSortedWl = sort(dMeanWl(:));
                        nnWl = numel(dSortedWl);
                        if nnWl == 0
                            continue;
                        end
                        idxLo = max(1, round(0.01 * nnWl)); %取 1% 分位数作为低端 lo，取 99% 分位数作为高端 hi
                        idxHi = max(1, round(0.99 * nnWl));
                        lo = dSortedWl(idxLo);
                        hi = dSortedWl(idxHi);
                        dRangeByWl(wi,:) = [max(lo*0.1, eps), hi*10];
                    end
                end
                if isempty(dRangeByWl)
                    dRangeByWl = params.fixedDRange;
                    wlVals = 1;
                    note = note + "[dRangeByWl auto fallback=single] ";
                end
            else
                if ~isempty(wl) && numel(wl) == nMeas
                    wlVals = unique(wl);
                    dRangeByWl = repmat(params.fixedDRange, numel(wlVals), 1);
                else
                    dRangeByWl = params.fixedDRange;
                    wlVals = 1;
                    note = note + "[dRangeByWl fallback=single] ";
                end
            end

            if ~isempty(dRangeByWl)
                fprintf('dRange used (by wl):\n');
                for wi = 1:numel(wlVals)
                    fprintf('  wl %d: [%.3g, %.3g]\n', wlVals(wi), dRangeByWl(wi,1), dRangeByWl(wi,2));
                end
            end

            %% ========== STEP 4: PRUNE CHANNELS (mlActAuto) ==========
            mlActAuto = [];
            usedPrune = false;
            nFailSNRMeas = NaN;
            nFailSDMeas = NaN;
            nFailSCIMeas = NaN;
            nFailSNRPair = NaN;
            nFailSDPair = NaN;
            nFailSCIPair = NaN;

            % Always use custom prune to apply per-wavelength dRange.
            note = note + "[PruneChannels custom per-wl] ";

            if usedPrune && ~iscell(mlActAuto)
                mlActAuto = {mlActAuto};
            end

            if usedPrune
                if isempty(mlActAuto) || isempty(mlActAuto{1}) || numel(mlActAuto{1}) ~= nMeas
                    warning('PruneChannels output size mismatch; using custom mask.');
                    note = note + "[PruneChannels size mismatch->custom] ";
                    usedPrune = false;
                end
            end

            if ~usedPrune
                % Custom prune for multi-wavelength (e.g., 3 WL)
                % hmrR_PruneChannels 只适用于2WL
                meanI = mean(d, 1, 'omitnan')';
                stdI  = std(d, 0, 1, 'omitnan')';
                snr   = meanI ./ (stdI + eps);

                if ~isempty(dRangeByWl) && ~isempty(wl) && numel(wl) == nMeas
                    inRange = false(nMeas,1);
                    for wi = 1:numel(wlVals)
                        idx = (wl == wlVals(wi));
                        inRange(idx) = (meanI(idx) >= dRangeByWl(wi,1)) & (meanI(idx) <= dRangeByWl(wi,2));
                    end
                else
                    inRange = (meanI >= dRangeByWl(1,1)) & (meanI <= dRangeByWl(1,2));
                end
                snrOK   = snr >= params.SNRthresh;
                sdOK    = true(nMeas,1);
                if ~all(isnan(dist))
                    sdOK = (dist >= params.SDrange(1)) & (dist <= params.SDrange(2));
                end
                measOK = inRange & snrOK & sdOK & sciOK;
                nFailSNRMeas = sum(~snrOK);
                nFailSDMeas = sum(~sdOK);
                if params.useSCI && sciComputed
                    nFailSCIMeas = sum(~sciOK);
                end
                if ~isempty(src) && ~isempty(det) && numel(src) == nMeas
                    [gPair, ~] = findgroups(src, det);
                    nPair = max(gPair);
                    if nPair > 0
                        failSNRPair = false(nPair,1);
                        failSDPair = false(nPair,1);
                        failSCIPair = false(nPair,1);
                        for gi = 1:nPair
                            idx = (gPair == gi);
                            failSNRPair(gi) = any(~snrOK(idx));
                            failSDPair(gi) = any(~sdOK(idx));
                            if params.useSCI && sciComputed
                                failSCIPair(gi) = any(~sciOK(idx));
                            end
                        end
                        nFailSNRPair = sum(failSNRPair);
                        nFailSDPair = sum(failSDPair);
                        if params.useSCI && sciComputed
                            nFailSCIPair = sum(failSCIPair);
                        end
                    end
                end

                if ~isempty(src) && ~isempty(det) && ~isempty(wl) && numel(src) == nMeas
                    mlActVec = zeros(nMeas,1);
                    [g, ~] = findgroups(src, det);
                    for gi = 1:max(g)
                        idx = (g == gi);
                        if numel(unique(wl(idx))) == nW && all(measOK(idx))
                            mlActVec(idx) = 1;
                        end
                    end
                else
                    mlActVec = measOK;
                end
                mlActAuto = {mlActVec};
                if ~isempty(nW) && nW ~= 2
                    note = note + "[PruneChannels custom nW!=2] ";
                else
                    note = note + "[PruneChannels custom] ";
                end
            end

            mlActVec = mlActAuto{1};
            nMeasTotal  = numel(mlActVec);
            nMeasActive = sum(mlActVec(:) ~= 0);
            activeMeasFrac = nMeasActive / max(nMeasTotal,1);

            fprintf('Active measurements: %d/%d (%.1f%%)\n', nMeasActive, nMeasTotal, 100*activeMeasFrac);

            %% ========== STEP 5: INTENSITY -> OD ==========
            dod = hmrR_Intensity2OD(intensity);

            %% ========== STEP 6: MOTION DETECTION ON OD (tIncAuto/tIncAutoCh) ==========
            tIncAuto = {ones(nT,1)};
            tIncAutoCh = {[]};
            flag_motionDetectFailed = false;
            flag_motionDetectCritical = false;
            try
                [tIncAuto, tIncAutoCh] = hmrR_MotionArtifactByChannel( ...
                    dod, probe, mlActMan_cell, mlActAuto, tIncMan_cell, ...
                    params.tMotion, params.tMask, params.STDEVthresh, params.AMPthresh);
            catch ME2
                warning('Motion detection failed; using all-ones tInc. Reason: %s', ME2.message);
                note = note + "[MotionDetect failed->skip] ";
                flag_motionDetectFailed = true;
                flag_motionDetectCritical = true;
            end

            if ~iscell(tIncAuto)
                tIncAuto = {tIncAuto};
            end
            if ~iscell(tIncAutoCh)
                tIncAutoCh = {tIncAutoCh};
            end

            tIncVec = tIncAuto{1};
            if isempty(tIncVec) || (~isnumeric(tIncVec) && ~islogical(tIncVec))
                warning('tIncAuto is invalid; resetting to all ones.');
                tIncVec = ones(nT,1);
                tIncAuto = {tIncVec};
                note = note + "[MotionDetect invalid->reset] ";
            end
            tIncVec = double(tIncVec(:));
            motionFrac = NaN;
            if ~flag_motionDetectFailed
                tIncCh = [];
                if ~isempty(tIncAutoCh) && ~isempty(tIncAutoCh{1})
                    tIncCh = tIncAutoCh{1};
                    if size(tIncCh,1) ~= nT && size(tIncCh,2) == nT
                        tIncCh = tIncCh';
                    end
                    if size(tIncCh,1) ~= nT
                        if size(tIncCh,1) > nT
                            tIncCh = tIncCh(1:nT, :);
                            note = note + "[MotionDetect tIncCh trimmed] ";
                        else
                            padN = nT - size(tIncCh,1);
                            tIncCh = [tIncCh; ones(padN, size(tIncCh,2))];
                            note = note + "[MotionDetect tIncCh padded] ";
                        end
                    end
                end
                if ~isempty(tIncCh) && size(tIncCh,1) == nT
                    tIncCh = double(tIncCh > 0);
                    tIncAutoCh{1} = tIncCh;
                    mlActVec = mlActAuto{1};
                    activeIdx = (mlActVec(:) ~= 0);
                    tIncChUse = tIncCh;
                    if any(activeIdx) && size(tIncChUse,2) == numel(activeIdx)
                        tIncChUse = tIncChUse(:, activeIdx);
                    end
                    motionFrac = 1 - mean(tIncChUse(:));
                else
                    if isempty(tIncCh)
                        flag_motionDetectCritical = true;
                        note = note + "[MotionDetect tIncCh empty] ";
                    else
                        flag_motionDetectCritical = true;
                        note = note + "[MotionDetect tIncCh dim mismatch] ";
                    end
                    motionFrac = 1 - mean(tIncVec);
                end
            end

            fprintf('Motion fraction (tInc==0): %.1f%%\n', 100*motionFrac);

            % Ensure mlActAuto is cell for Homer3 motion correction
            if ~iscell(mlActAuto)
                mlActAuto = {mlActAuto};
            end

            %% ========== STEP 7: MOTION CORRECTION ==========
            % 7a) TDDR (Temporal Derivative Distribution Repair)
            didTDDR = false;
            if params.doTDDR && tddrAvailable
                try
                    dt = median(diff(t));
                    if ~isfinite(dt) || dt <= 0
                        error('TDDRInvalidTime:NonPositiveDT', 'Invalid time vector for TDDR.');
                    end
                    fs = 1 / dt;
                    dodData = dod.dataTimeSeries;
                    if ~ismatrix(dodData)
                        error('TDDRInvalidData:DataTimeSeries', 'dod.dataTimeSeries has invalid shape.');
                    end
                    validCols = all(isfinite(dodData), 1);
                    if any(validCols)
                        dodData(:, validCols) = TDDR(dodData(:, validCols), fs);
                    end
                    dod.dataTimeSeries = dodData;
                    didTDDR = true;
                catch MEtd
                    error('TDDRFailed:MotionCorrectTDDR', ...
                        'TDDR motion correction failed: %s', MEtd.message);
                end
            else
                note = note + "[TDDR disabled] ";
            end

            % 7b) Wavelet (recommended, but may fail if wavelet toolbox missing) 
            % 这一步时间比较长
            didWavelet = false;
            if params.doWavelet && waveletAvailable
                try
                    dod = hmrR_MotionCorrectWavelet(dod, mlActMan_cell, mlActAuto, params.iqrWavelet, 1);
                    didWavelet = true;
                catch MEwv
                    error('WaveletFailed:MotionCorrectWavelet', ...
                        'Wavelet motion correction failed: %s', MEwv.message);
                end
            else
                note = note + "[Wavelet disabled] ";
            end

            %% ========== STEP 8: BANDPASS FILTER (OD) ==========
            try
                dod = hmrR_BandpassFilt(dod, params.hpf, params.lpf);
            catch MEbp
                error('Bandpass filter failed: %s', MEbp.message);
            end

            %% ========== STEP 9: OD -> Hb (OD2Conc) WITH 3 WAVELENGTHS ==========
            % PPF must have length = nW. We replicate scalar to vector.
            ppfUsed = repmat(params.ppfScalar, 1, nW);

            % Some datasets might report nW but still require at least 2; you have 3 -> OK.
            dc = hmrR_OD2Conc(dod, probe, ppfUsed);

            %% ========== STEP 10: QC FLAGS ==========
            flag_lowActive  = (activeMeasFrac < qcRule.minActiveMeasFrac);
            flag_highMotion = (motionFrac > qcRule.maxMotionFrac);

            %% ========== STEP 11: SAVE MAT ==========
            outMat = fullfile(subPath, sprintf('%s_task-audi_fnirs_preproc.mat', subName));
            qc = struct();
            qc.sub = subName;
            qc.run = runBase;
            qc.inFile = inFile;
            qc.durationSec = durationSec;
            qc.nW = nW;
            qc.nMeasTotal = nMeasTotal;
            qc.nMeasActive = nMeasActive;
            qc.activeMeasFrac = activeMeasFrac;
            qc.motionFrac = motionFrac;
            qc.nonPositiveMeasCount = nonPositiveMeasCount;
            qc.nonPositiveMeasFrac = nonPositiveMeasFrac;
            qc.nFailSNRMeas = nFailSNRMeas;
            qc.nFailSDMeas = nFailSDMeas;
            qc.nFailSCIMeas = nFailSCIMeas;
            qc.nFailSNRPair = nFailSNRPair;
            qc.nFailSDPair = nFailSDPair;
            qc.nFailSCIPair = nFailSCIPair;
            qc.flag_nonPositiveCritical = flag_nonPositiveCritical;
            qc.flag_nonPositiveSevere = flag_nonPositiveSevere;
            qc.flag_motionDetectFailed = flag_motionDetectFailed;
            qc.flag_motionDetectCritical = flag_motionDetectCritical;
            qc.didTDDR = didTDDR;
            qc.didWavelet = didWavelet;
            qc.sciPairs = sciPairs;
            qc.sciMed = sciMed;
            qc.sciGoodPct = sciGoodPct;
            qc.sciThr = params.SCIthr;
            qc.sciWlPair = params.SCIwlPair;
            qc.dRangeByWl = dRangeByWl;
            qc.dRangeWlIndex = wlVals;
            qc.ppfUsed = ppfUsed;
            qc.addedOffset = addedOffset;
            qc.flag_lowActive = flag_lowActive;
            qc.flag_highMotion = flag_highMotion;
            qc.note = char(note);


            stim = [];

            stim = snirfObj.stim; 
            % Save the key intermediate masks too (useful for GLM later)
            save(outMat, 'dod', 'dc', 'qc', 'params', 'mlActAuto', 'tIncAuto', 'tIncAutoCh', 'stim', '-v7.3');

            %% ========== STEP 12: QC FIGURES ==========
            % (Plotting removed by request.)

            %% ========== STEP 13: ADD QC TABLE ROW ==========
            dRangeByWlStr = '';
            if ~isempty(dRangeByWl)
                for wi = 1:numel(wlVals)
                    dRangeByWlStr = [dRangeByWlStr, sprintf('wl%d:[%.3g %.3g] ', wlVals(wi), dRangeByWl(wi,1), dRangeByWl(wi,2))]; %#ok<AGROW>
                end
                dRangeByWlStr = strtrim(dRangeByWlStr);
            end

            sciMedStr = '';
            sciGoodPctStr = '';
            if ~isempty(sciPairs) && ~isempty(sciMed)
                nPairs = size(sciPairs,1);
                for pi = 1:nPairs
                    if pi <= numel(sciMed) && isfinite(sciMed(pi))
                        sciMedStr = [sciMedStr, sprintf('s%d-d%d:%.3g ', sciPairs(pi,1), sciPairs(pi,2), sciMed(pi))]; %#ok<AGROW>
                    end
                    if pi <= numel(sciGoodPct) && isfinite(sciGoodPct(pi))
                        sciGoodPctStr = [sciGoodPctStr, sprintf('s%d-d%d:%.3g ', sciPairs(pi,1), sciPairs(pi,2), sciGoodPct(pi))]; %#ok<AGROW>
                    end
                end
                sciMedStr = strtrim(sciMedStr);
                sciGoodPctStr = strtrim(sciGoodPctStr);
            end

            sciWlPairStr = '';
            if ~isempty(params.SCIwlPair)
                sciWlPairStr = sprintf('wl%d-wl%d', params.SCIwlPair(1), params.SCIwlPair(2));
            end

            qcRow = {subName, runBase, inFile, outMat, nW, durationSec, nMeasTotal, nMeasActive, ...
                     activeMeasFrac, motionFrac, nonPositiveMeasFrac, nonPositiveMeasCount, ...
                     nFailSNRMeas, nFailSDMeas, nFailSCIMeas, nFailSNRPair, nFailSDPair, nFailSCIPair, ...
                     flag_nonPositiveCritical, flag_nonPositiveSevere, flag_motionDetectFailed, ...
                     flag_motionDetectCritical, didTDDR, didWavelet, sciMedStr, sciGoodPctStr, ...
                     params.SCIthr, sciWlPairStr, dRangeByWlStr, ...
                     strtrim(sprintf('%.2g ', ppfUsed)), addedOffset, flag_lowActive, ...
                     flag_highMotion, char(note)};
            qcData = [qcData; qcRow]; %#ok<AGROW>

            fprintf('[DONE] %s | nW=%d | active=%.1f%% | motion=%.1f%%\n', ...
                runTag, nW, 100*activeMeasFrac, 100*motionFrac);

        catch ME
            if strcmp(ME.identifier, 'WaveletFailed:MotionCorrectWavelet') || ...
               strcmp(ME.identifier, 'TDDRFailed:MotionCorrectTDDR')
                rethrow(ME);
            end
            % If any run fails, we still write a QC row with the error message.
            warning('[ERROR] %s failed: %s', runTag, ME.message);

            qcRow = {subName, runBase, inFile, '', NaN, NaN, NaN, NaN, ...
                     NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                     NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                     '', '', NaN, '', '', '', NaN, NaN, NaN, ...
                     ['ERROR: ' ME.message]};
            qcData = [qcData; qcRow]; %#ok<AGROW>
        end

        if debug.processOnlyFirstFile
            break;
        end
    end
end

%% ===================== [5] WRITE QC SUMMARY CSV =====================
if ~isempty(qcData)
    T = cell2table(qcData, 'VariableNames', qcHeader);
    outCsv = fullfile(outDir, 'QC_summary.csv');
    writetable(T, outCsv);
    fprintf('\nQC summary saved: %s\n', outCsv);
else
    fprintf('\nNo runs processed. Check rootDir and SNIRF files.\n');
end
