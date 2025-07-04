# EEG Sleep Spindle Analysis in COVID-19 Survivors

This academic project investigates sleep spindles in a COVID-19 ICU survivor and a control subject using high-density EEG data. The analysis compares slow and fast spindles in terms of spectral properties, topographic distribution, and inverse localization using EEGLAB and Brainstorm.

## Project Description
This project was carried out for the Neurorobotics and Neurorehabilitation course at the University of Padova. The aim was to evaluate changes in neurophysiological sleep markers in post-COVID individuals and assess potential evidence of stress-induced cognitive impairment.

## Aims of the Study
- Identify and classify slow (9–12 Hz) and fast (12–16 Hz) spindles from nap EEG recordings.
- Compute power spectra and discard out-of-range spindles.
- Generate topographic maps of spindle activity using EEGLAB.
- Localize spindle sources at 0 ms using sLORETA in Brainstorm.

## Technologies & Tools
- Language: MATLAB
- Tools & Toolboxes:
  - EEGLAB (for topographies)
  - Brainstorm (for source localization)
  - Signal Processing Toolbox (butter, filtfilt, pwelch)
  - Custom scripts for spindle extraction and spectrum filtering

The dataset used is too large to upload here but is available [here](https://drive.google.com/drive/folders/1t0kmrh_W7CI5zuTup-H1F2yx1g7dt-E5?usp=share_link).

## Folder Structure
sleep-spindle-analysis-covid/
│
├── main_script.m                # MATLAB script implementing the full analysis pipeline
├── solution_report.pdf          # Final project report with results and discussions
├── assignment.pdf               # Original project assignment
└── README.md


## Main Steps of the Analysis

### 1. Data Preparation
- Load nap EEG recordings (204 channels, 250 Hz).
- Generate time vectors for both subjects.

### 2. Filtering
- Apply band-pass filters to isolate slow (9–12 Hz) and fast (12–16 Hz) spindles.

### 3. Spindle Extraction
- Use fixed 500 ms windows starting from provided spindle timepoints.
- Build 3D matrices of spindle signals.

### 4. Channel Averaging
- Average EEG channels per spindle to obtain 1D signals for spectral analysis.

### 5–6. Spectral Analysis & Selection
- Compute power spectra using Welch's method.
- Retain only those with peak frequencies within defined ranges.

### 7. Averaging of Selected Spindles
- Average across selected spindle trials and channels.

### 8. Topographic Mapping
- Plot topographies using EEGLAB to visualize spatial distribution.

### 9. Source Localization
- Import averaged data into Brainstorm.
- Use sLORETA to localize sources at 0 ms in inverse space.

## Key Findings
- Slow spindles showed higher spectral amplitude than fast ones.
- ICU subject’s spindles were more anterior and less powerful than those of the control.
- Activity localization revealed neurophysiological differences that may reflect post-traumatic stress and cognitive changes following COVID-19 ICU admission.
