<!-- Copilot / AI agent instructions for the P3 EEG ERP Processing repo -->
# P3 EEG ERP Processing — Copilot instructions

Purpose: give an AI agent the minimum, actionable context to reason about and modify this MATLAB-based ERP pipeline safely.

- Project type: MATLAB scripts that rely on EEGLAB + ERPLAB toolboxes. Primary processing is implemented as numbered scripts: Script1_... through Script8_... (per-subject), plus downstream measurement and grand-average scripts in `ERP_Measurements` and `Grand_Average_ERPs`.

- Big picture / pipeline (order matters):
  - `Script1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt.m` — import .set, shift event codes, downsample, rereference, create bipolars, add locs, HP filter.
  - `Script2_ICA_Prep.m` — remove breaks/noisy segments and prepare files for ICA (reads `ICA_Prep_Values_P3.xls`).
  - `Script3_Run_ICA.m` — (commented by default) compute ICA weights if needed; note: recomputing ICA changes component ordering.
  - `Script4_Remove_ICA_Components.m` — load ICA weights and remove components listed in `ICA_Components_P3.xlsx`.
  - `Script5_Elist_Bin_Epoch.m` — make event list, bin with `BDF_P3.txt`, epoch and baseline-correct.
  - `Script6_Artifact_Rejection.m` — interpolate channels, run multiple artifact detectors (reads several per-subject Excel sheets), save cleaned sets.
  - `Script7_Average_ERPs.m` — average ERPs, compute percent-rejected, create diff waves and low-pass filtered versions.
  - `Script8_Plot_Individual_Subject_ERPs.m` — generate per-subject PDFs of plots.
  - `ERP_Measurements/Script12_Measure_ERPs.m` — extract measurement values (mean, peak, latencies) using lists created by Script7.
  - `Grand_Average_ERPs/Script9_Grand_Average_ERPs.m` and plotting scripts (`Script10`, `Script11`) — build and plot grand averages.

- Data & naming conventions (important):
  - Subject folders are top-level numeric names (e.g., `1`, `2`, ...). Many scripts build filenames like `SUBID_P3_*` in each subject folder.
  - Intermediate .set and .erp filenames follow strict suffix patterns (e.g., `_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip.set`, `_P3_erp_ar_lpfilt.erp`) — preserve naming when producing downstream inputs.
  - Per-subject configuration is provided in Excel files and plain text BDF files at repo root: `ICA_Components_P3.xlsx`, `ICA_Prep_Values_P3.xls`, `Interpolate_Channels_P3.xls`, `AR_Parameters_for_*.xls`, `BDF_P3.txt`, `P3_Diff_Wave.txt`, `Rereference_Add_Uncorrected_Bipolars_P3.txt`, `Add_Corrected_Bipolars_P3.txt`.

- Tooling & workflows:
  - Execution environment: run interactively in MATLAB (or Octave if compatible) with EEGLAB and ERPLAB on the path. Scripts call `eeglab` and many `pop_*` / `erp_*` functions.
  - Full pipeline: open MATLAB, ensure current folder is `EEG_ERP_Processing`, then run the numbered scripts in order (1 → 8). Do not run Script3 by default—it's intentionally commented out in upstream distribution.
  - To re-run only later stages for a subset of subjects: edit the `SUB` variable at top of each script (it lists subject IDs) or run the script for a single `SUB{i}` inside MATLAB.

- Project-specific patterns & safe-edit rules for AI agents:
  - These are procedural scripts (not functions). Avoid converting scripts into functions without explicit user approval — users rely on workspace variables and file-system side effects.
  - Do not change filenames or output locations unless you update every downstream script that depends on them. Prefer adding optional parameters or wrapper scripts.
  - Many per-subject thresholds are stored in Excel sheets. For behavioral/parameter changes, prefer editing the Excel inputs rather than hard-coding new constants into scripts.
  - `Script3_Run_ICA.m` contains a clear warning: recomputing ICA will invalidate `ICA_Components_P3.xlsx`. If an agent recomputes ICA, it must also update component-selection records and inform the user.

- Integration points & external dependencies:
  - EEGLAB and ERPLAB toolboxes (MATLAB). Some code references `binica` but falls back to `runica` — be conservative when changing ICA calls.
  - Excel files are read with `xlsread` — keep the exact expected layout (subject ID in first column, parameters in following columns).
  - Outputs include pdfs created with `save2pdf` and .csv/.txt measurement exports; preserve these for reproducibility.

- Guidance for code edits and PRs:
  - Small, localized changes are preferred (fix a path, correct a filename, update one parameter file). When making pipeline/format changes, include a migration plan that updates all affected scripts.
  - Add tests only after confirming how the user runs code (no test harness exists). Ask before creating CI or converting to functions.

Examples (explicit):
- To run the whole pipeline interactively: open MATLAB, set current folder to `EEG_ERP_Processing`, then run each script in order: `Script1_...` → ... → `Script8_...`.
- To change which subjects are processed, edit the `SUB` cell array at the top of any `Script*.m`.

If anything here is unclear or you want the agent to follow stricter rules (for example: always open a draft PR, never recompute ICA, or prefer adding wrapper functions), tell me which rule to enforce and I will update this file.
