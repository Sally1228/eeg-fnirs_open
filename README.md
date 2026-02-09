# EEG-fNIRS Analysis Pipeline for Cerebral Small Vessel Disease

**English** | [中文](README_CN.md)

A comprehensive MATLAB-based pipeline for simultaneous EEG-fNIRS data processing, feature extraction, and neurovascular coupling analysis in cerebral small vessel disease (CSVD) subtypes.

## Overview

This project implements a complete multimodal neuroimaging analysis framework to investigate brain functional differences across three CSVD subtypes — sporadic SVD, CADASIL, and CAA — compared with healthy controls (HC). The pipeline covers three cognitive-motor paradigms and integrates EEG electrophysiology with fNIRS hemodynamics to characterize neurovascular coupling impairments.

### Task Paradigms

| Task | Paradigm | EEG Features | fNIRS ROIs |
|------|----------|-------------|-----------|
| **Auditory Oddball** | Standard/deviant tone detection | P300 ERP, Delta/Theta power | PFC, Parietal, Temporal |
| **Finger Tapping** | Left/right hand motor execution | Mu/Beta ERD, lateralization index | Bilateral central regions |
| **Semantic Categorization** | Visual word categorization | Alpha/Beta relative power, ERD | Frontal, Temporal |

## Repository Structure

```
├── eeg/                                # EEG data & processing
│   ├── audi/                           # Auditory oddball task
│   │   ├── audi_EEG_ERP_Processing/
│   │   │   └── scripts/               # STEP0–STEP14 MATLAB scripts
│   │   └── sub-*/                      # Per-subject raw data (.edf)
│   ├── motor/                          # Finger tapping task
│   │   ├── motor_EEG_Processing/
│   │   │   └── scripts/               # STEP0–STEP7 MATLAB scripts
│   │   └── sub-*/
│   └── cate/                           # Semantic categorization task
│       ├── cate_EEG_Processing/
│       │   └── scripts/               # STEP0–STEP7 MATLAB scripts
│       └── sub-*/
│
├── fnirs/                              # fNIRS data & processing
│   ├── audi/
│   │   ├── audi_fnirs_processing/     # STEP1–STEP6 MATLAB scripts
│   │   └── sub-*/                      # Per-subject data (.snirf)
│   ├── motor/
│   │   ├── motor_fnirs_processing/
│   │   └── sub-*/
│   ├── cate/
│   │   ├── cate_fnirs_processing/
│   │   └── sub-*/
│   └── nirs_to_snirf.m               # Format conversion utility
│
└── eeg_fnirs_processing/              # EEG-fNIRS coupling analysis
    ├── audi_eeg_fnirs_processing/     # STEP1–STEP4 (full coupling pipeline)
    ├── motor_eeg_fnirs_processing/    # STEP1–STEP2 (block-level + EEG-informed)
    └── cate_eeg_fnirs_processing/     # STEP1–STEP2 (block-level + EEG-informed)
```

## Processing Pipelines

### 1. EEG Processing

Sequential MATLAB scripts following the EEGLAB/ERPLAB framework:

| Step | Script | Description |
|------|--------|-------------|
| 0 | `STEP0_marker.m` | Extract event markers from trigger channel; correct LCD monitor delay |
| 1 | `STEP1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt.m` | Import raw EEG (1000 Hz) → downsample to 256 Hz → re-reference to averaged mastoids (A1/A2) → 0.1 Hz high-pass filter |
| 2 | `STEP2_ICA_Prep.m` | Prepare data for ICA (channel interpolation, epoch rejection) |
| 3 | `STEP3_Run_ICA.m` | Extended Infomax ICA decomposition (60-channel) with ICLabel auto-classification |
| 4 | `STEP4_Remove_ICA_Components.m` | Remove ocular/muscular artifact components |
| 5 | `STEP5_Elist_Bin_Epoch.m` | Create event list, assign bins, epoch extraction (−200 to 2000 ms) |
| 6 | `STEP6_Artifact_Rejection.m` | Artifact detection: simple voltage threshold + moving-window peak-to-peak |
| 7+ | Task-specific | **Oddball**: `STEP7_Average_ERPs.m` → ERP averaging and measurement; **Motor**: `STEP7_ERD.m` → time-frequency ERD analysis (wavelet); **Semantic**: `STEP7_power.m` → spectral power extraction |

### 2. fNIRS Processing

Sequential MATLAB scripts based on the Homer3 framework:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `STEP1_preprocessing.m` | Load SNIRF → convert to OD → channel quality assessment (SNR ≥ 2, SCI ≥ 0.6) → motion artifact correction (TDDR + wavelet, IQR = 1.5) → bandpass filter (0.01–0.2 Hz, Butterworth) → modified Beer-Lambert Law (PPF = 6.0) → output ΔHbO/ΔHbR |
| 2 | `STEP2_qc.m` | Quality control visualization: raw intensity → OD → corrected OD → HbO/HbR time courses; generate `QC_summary.csv` |
| 3 | `STEP3_BlockAverage.m` | Block averaging: identify task blocks → baseline correction (−5 to 0 s) → extract block-averaged hemodynamic responses |
| 4 | `STEP4_GLM.m` | Subject-level GLM: stimulus × gamma HRF (k=6, θ=1, 30 s duration) convolution → design matrix with DCT drift regressors (256 s cutoff) → OLS/AR-IRLS fitting → output β coefficients and t-statistics |
| 5 | `STEP5_group_block.m` | Group-level block analysis with one-sample t-test and FDR correction |
| 6 | `STEP6_group_GLM.m` | Group-level GLM analysis |

### 3. EEG-fNIRS Coupling Analysis

Multi-level neurovascular coupling analysis combining both modalities:

| Step | Script | Method |
|------|--------|--------|
| 1 | `STEP1_correlation.m` | **Pearson correlation** between EEG features (P300 amplitude, band power, ERD) and fNIRS hemodynamics (ΔHbO/ΔHbR) per ROI. Oddball: trial-level; Motor/Semantic: block-level |
| 2 | `STEP2_EEG_informed_fnirs.m` | **EEG-informed GLM**: EEG features (z-scored) modulate HRF → parametric regressor in fNIRS GLM → quantify neural activity's contribution to hemodynamic variance |
| 3 | `STEP3_timelag_correlation_and_coherence.m` | **Cross-correlation & coherence**: Hilbert envelope of EEG band power → interpolate to fNIRS sampling rate → xcorr for optimal lag time → magnitude-squared coherence (Oddball task only) |
| 4 | `STEP4_dynamic_coupling.m` | **Dynamic coupling models**: ARX system identification (na=nb=2, delay 2.5 s) → state-space modeling → adaptive Kalman filtering (forgetting factor 0.98) for time-varying coupling estimation (Oddball task only) |

## Key Parameters

### EEG

| Parameter | Value |
|-----------|-------|
| Sampling rate (raw → processed) | 1000 Hz → 256 Hz |
| Reference | Averaged mastoids (A1/A2) |
| High-pass filter | 0.1 Hz |
| Low-pass filter (ERP) | 20 Hz |
| Epoch window | −200 to 2000 ms |
| Baseline correction | −200 to 0 ms |
| ICA algorithm | Extended Infomax |
| P300 time window | 300–500 ms at Cz/Pz |

### fNIRS

| Parameter | Value |
|-----------|-------|
| Channel quality | SNR ≥ 2, SCI ≥ 0.6 |
| Motion correction | TDDR + wavelet denoising (IQR = 1.5) |
| Bandpass filter | 0.01–0.2 Hz (Butterworth) |
| Beer-Lambert PPF | 6.0 |
| HRF model | Gamma function (k = 6, θ = 1, 30 s) |
| GLM drift removal | DCT high-pass (cutoff = 256 s) |
| Block baseline | −5 to 0 s |

## Data Format

| Modality | Raw Format | Processed Format |
|----------|-----------|-----------------|
| EEG | `.edf` (European Data Format) | `.set` (EEGLAB) |
| fNIRS | `.snirf` (Shared Near-Infrared Spectroscopy Format) | `.mat` (MATLAB) |
| Behavioral | `.txt` (E-Prime export) | `.csv` |
| Results | — | `.csv`, `.mat`, `.xlsx` |

## Dependencies

- **MATLAB** R2020b or later
- **EEGLAB** (v2021+) with plugins:
  - ERPLAB
  - ICLabel
  - clean_rawdata
- **Homer3** (fNIRS processing toolbox)
- **Signal Processing Toolbox** (MATLAB)
- **Statistics and Machine Learning Toolbox** (MATLAB)

## Usage

Each processing module follows a sequential STEP-based workflow. Scripts auto-discover subject directories matching the `sub-*` pattern.

```matlab
% Example: Run EEG preprocessing for the auditory oddball task
cd eeg/audi/audi_EEG_ERP_Processing/scripts/
STEP0_marker
STEP1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt
STEP2_ICA_Prep
STEP3_Run_ICA
STEP4_Remove_ICA_Components   % requires manual inspection
STEP5_Elist_Bin_Epoch
STEP6_Artifact_Rejection
STEP7_Average_ERPs

% Example: Run fNIRS preprocessing
cd fnirs/audi/audi_fnirs_processing/
STEP1_preprocessing
STEP2_qc
STEP3_BlockAvarage
STEP4_GLM

% Example: Run EEG-fNIRS coupling analysis
cd eeg_fnirs_processing/audi_eeg_fnirs_processing/
STEP1_correlation
STEP2_EEG_informed_fnirs
STEP3_timelag_correlation_and_coherence
STEP4_dynamic_coupling
```

> **Note:** Steps must be executed in order. `STEP4_Remove_ICA_Components` in EEG processing requires manual visual inspection of ICA components.

## Analysis Flowchart

```
Raw EEG (.edf)                     Raw fNIRS (.snirf)
      │                                   │
      ▼                                   ▼
┌─────────────────┐             ┌──────────────────────┐
│  EEG Processing │             │  fNIRS Processing    │
│  STEP 0–6       │             │  STEP 1–4            │
│  ·Resampling    │             │  ·Quality control    │
│  ·Re-reference  │             │  ·Motion correction  │
│  ·ICA           │             │  ·Bandpass filtering │
│  ·Artifact rej. │             │  ·MBLL conversion    │
└────────┬────────┘             └──────────┬───────────┘
         │                                 │
         ▼                                 ▼
┌─────────────────┐             ┌──────────────────────┐
│ Feature Extract  │             │ Feature Extraction   │
│ ·P300 (Oddball) │             │ ·Block avg ΔHbO/HbR │
│ ·ERD (Motor)    │             │ ·GLM β coefficients  │
│ ·Power (Cate)   │             │ ·Peak time           │
└────────┬────────┘             └──────────┬───────────┘
         │                                 │
         └──────────┬──────────────────────┘
                    ▼
    ┌───────────────────────────────────┐
    │  EEG-fNIRS Coupling Analysis      │
    │  ·Block/trial-level correlation   │
    │  ·EEG-informed GLM               │
    │  ·Cross-correlation & coherence   │
    │  ·ARX / Kalman dynamic modeling   │
    └───────────────┬───────────────────┘
                    ▼
    ┌───────────────────────────────────┐
    │  Group-level Statistics           │
    │  ·Welch t-test (vs HC)           │
    │  ·Benjamini-Hochberg FDR         │
    └───────────────────────────────────┘
```

## Citation

If you use this pipeline in your research, please cite:

> [Your publication information here]

## License

[License information here]
