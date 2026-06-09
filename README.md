# uSMAART — ultra-Sensitive Measurement of Aggregate Activity in Restricted cell-Types

MATLAB toolkit for preprocessing, unmixing, and analyzing fiber photometry recordings acquired with the uSMAART instrument. uSMAART is a dual-wavelength laser fiber-photometry system designed to detect high-frequency voltage dynamics from genetically identified cell types using GEVIs (genetically encoded fluorescent voltage indicators). The instrument supports up to **2 sensors** recording from up to **2 brain areas** simultaneously.

> **Reference:** Haziza S. et al. *Imaging high-frequency voltage dynamics in multiple neuron classes of behaving mammals.* Cell, 2025. https://doi.org/10.1016/j.cell.2025.06.028

---

## Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Repository Structure](#repository-structure)
5. [Typical Processing Pipeline](#typical-processing-pipeline)
6. [Function Reference](#function-reference)
7. [Demo Files](#demo-files)
8. [Known Limitations & Notes](#known-limitations--notes)
9. [License](#license)

---

## Requirements

- **MATLAB R2019b or later** (tested up to R2024a)
- **MATLAB Toolboxes:**
  - Signal Processing Toolbox (required — filtering, spectral analysis)
  - Statistics and Machine Learning Toolbox (required — `robustfit`, `fitlm`)
  - Wavelet Toolbox (optional — `powerERP` time-frequency analysis)
- **Third-party:**
  - [FastICA](https://research.ics.aalto.fi/ica/fastica/) — required only for `umxICA` (ICA-based unmixing). Download and add to your MATLAB path before using this function.

---

## Installation

```matlab
% 1. Clone the repository
git clone https://github.com/sihaziza/uSMAART_public.git
cd uSMAART_public

% 2. Add all subfolders to the MATLAB path (run once per session, or save permanently)
addpath(genpath(pwd))
savepath   % optional: saves for future sessions
```

---

## Quick Start

Open one of the four demo scripts in MATLAB. Each demo represents a different acquisition configuration:

| Demo file | Configuration | Description |
|---|---|---|
| `demo_1S1A.m` | 1 sensor, 1 area | Single GEVI, single recording site |
| `demo_2S1A.m` | 2 sensors, 1 area | Dual GEVI, single site (e.g., voltage + calcium) |
| `demo_2S2A.m` | 2 sensors, 2 areas | Dual GEVI, dual site |
| `demo_imaging.m` | Wide-field imaging | GEVI imaging from cortical window |

**Before running:** edit the `mainPath` / `mousePath` variables at the top of each demo to point to your data.

```matlab
% Minimal example: single sensor, single area
mainPath = 'path/to/your/data';
folder   = 'your_experiment_folder';
demo_1S1A
```

---

## Repository Structure

```
uSMAART_public/
│
├── demo_1S1A.m             — Demo: 1 sensor, 1 brain area
├── demo_2S1A.m             — Demo: 2 sensors, 1 brain area
├── demo_2S2A.m             — Demo: 2 sensors, 2 brain areas
├── demo_imaging.m          — Demo: wide-field imaging
│
├── preprocessing/
│   ├── loading/
│   │   ├── loaduSMAART2mat.m            — Load MFLI .mat files (standard)
│   │   └── loaduSMAART2mat_2CT2BR.m     — Load MFLI .mat files (2-color, 2-area)
│   │
│   ├── conditioning/
│   │   ├── bpFilter1D.m                 — Butterworth bandpass/highpass/lowpass (1D/2D/3D)
│   │   ├── sh_NotchFilter.m             — Alias for notchFilter (backwards compat.)
│   │   ├── runPhotoBleachingRemoval.m   — Photobleaching / drift correction
│   │   └── sh_zscore.m                  — Z-score normalization with range option
│   │
│   ├── modalityAlignment/
│   │   ├── alignMultimodalClockDrift.m  — Clock-drift correction via TTL pulses
│   │   ├── alignSyncTTL.m               — TTL cross-correlation alignment
│   │   ├── alignRecording.m             — Orchestrator: align ePhys + oPhys
│   │   ├── parseEPhys.m                 — Parse multi-epoch ePhys recordings
│   │   ├── parsingAligningMultimodalData.m — Multi-format data parser
│   │   └── read_Intan_RHD2000_file.m    — Intan RHD2000 binary file reader
│   │
│   ├── behaviorVideo/
│   │   ├── preDLCconditioning.m         — Pre-process AVI for DeepLabCut
│   │   ├── getMouseSpeedDLC.m           — Extract speed from DLC .h5 output
│   │   ├── getPixelCalibration.m        — Pixel-to-cm calibration
│   │   ├── estimateBehavioralState.m    — Wrapper: speed → behavioral state
│   │   ├── findRestRunTransitions.m     — Detect rest↔run transitions
│   │   └── getRestRunBootIndex.m        — Bootstrap indices for rest/run epochs
│   │
│   └── unmixing/
│       ├── unmixing1D.m                 — Router: select unmixing method
│       ├── umxCONV.m                    — Convolutional unmixing (Wiener filter)
│       ├── umxHDM.m                     — Heartbeat-driven minimization
│       ├── umxPCA.m                     — PCA-based unmixing
│       ├── umxICA.m                     — ICA-based unmixing (requires FastICA)
│       ├── umxRLR.m                     — Robust linear regression unmixing
│       ├── FindHBpeak.m                 — Automatic heartbeat frequency detection
│       ├── deconvolution.m              — Tikhonov-regularized deconvolution
│       └── estimateFilter.m             — Sliding-window filter estimation
│
├── analysis/
│   ├── Analysis_dualModality_OPhysEPhys.m  — Full oPhys + ePhys analysis pipeline
│   ├── Analysis_imaging.m                  — Wide-field imaging analysis pipeline
│   │
│   ├── eventRelatedAnalysis/
│   │   ├── getStimEpoch.m     — Extract stimulus-triggered epochs (core function)
│   │   ├── rasterERP.m        — Raster plot + event-related average
│   │   └── powerERP.m         — Time-frequency (CWT) event-related power
│   │
│   └── functions/
│       ├── notchFilter.m      — Notch filter (canonical implementation)
│       ├── findRipples.m      — Sharp-wave ripple detection
│       ├── findSpindles.m     — Sleep spindle detection
│       ├── cleanSWR.m         — Interactive ripple validation GUI
│       ├── coherencePair.m    — Signal pair coherence (mscohere)
│       ├── plotPSD.m          — Power spectral density visualization
│       └── plotErrorBar1.m    — Mean ± error band plot (SEM, std, CI95, CI99)
│
└── utilities/
    ├── getOptions.m           — Name-value pair argument parser
    ├── getTime.m              — Time vector from data and sampling rate
    ├── geviColor.m            — Color palette for GEVI types
    ├── istensor.m             — 3D array check
    ├── bpFilter2D.m           — 2D spatial Gaussian bandpass filter
    ├── convertBWtoF.m         — Bandwidth edges → frequency indices
    ├── getPointProjection.m   — Mean projection across spatial dimensions
    ├── getUnitsROI.m          — Extract ROI traces from H5 imaging data
    ├── get_circular_mask.m    — Interactive circular mask for GRIN lens imaging
    └── autoCropImage.m        — Intensity-threshold image cropping
```

---

## Typical Processing Pipeline

### Fiber photometry (1S1A / 2S1A / 2S2A)

```
1. Load raw MFLI data
   loaduSMAART2mat()

2. Signal conditioning
   sh_NotchFilter()            % remove power-line interference (e.g. 60 Hz or 297.6 Hz)
   bpFilter1D()                % bandpass to frequency range of interest
   runPhotoBleachingRemoval()  % correct slow drift / photobleaching

3. (If multi-modal) Align ePhys and oPhys clocks
   alignMultimodalClockDrift() + alignSyncTTL()

4. Hemodynamic unmixing  — choose one method:
   umxCONV()   % recommended default (convolutional, Wiener filter)
   umxRLR()    % robust linear regression (simpler, fast)
   umxPCA()    % PCA-based
   umxHDM()    % heartbeat-driven minimization
   umxICA()    % ICA (requires FastICA)

5. Event-triggered analysis
   rasterERP() / getStimEpoch() / powerERP()

6. Visualization
   plotPSD() / plotErrorBar1() / coherencePair()
```

### Wide-field imaging

```
1. Load H5 data, draw ROI
   getUnitsROI()

2. Load behavior
   getMouseSpeedDLC() → estimateBehavioralState()

3. Detrending
   runPhotoBleachingRemoval()

4. Unmixing
   umxCONV() with a range of epoch values; select best by PSD inspection

5. Event-triggered analysis + visualization
   rasterERP() / plotErrorBar1() / plotPSD()
```

---

## Function Reference

### Key parameters shared across most functions

| Parameter | Default | Description |
|---|---|---|
| `samplingRate` | 1000 | Acquisition sampling rate (Hz) |
| `VerboseMessage` | `true` | Print progress to command window |
| `VerboseFigure` | `true` | Generate diagnostic figures |

### Unmixing methods comparison

| Function | Method | Best used when |
|---|---|---|
| `umxCONV` | Wiener deconvolution | Default; handles most cases well |
| `umxRLR` | Robust linear regression | Simple, fast, low HB contamination |
| `umxPCA` | PCA | Strong linear mixing between channels |
| `umxHDM` | Heartbeat minimization | Strong cardiac artifact, clear HB peak |
| `umxICA` | FastICA | Complex mixing; requires FastICA installed |

### `getStimEpoch` — stimulus-triggered epoch extraction

```matlab
output = getStimEpoch(Data, TTL, Fs, 'baselinePrePost', 1, ...
                      'getShuffle', true, 'preNormalize', 'median')
```

`preNormalize` options: `'median'` (default), `'mean'`, `'zscore'`, `'none'`

Output fields: `arrayRaw`, `arrayShuffle`, `stimBand`, `indexTTLraw`, `indexTTLshuffle`

### `plotErrorBar1` — mean ± error visualization

```matlab
plotErrorBar1(data, 'x_axis', t, 'error', 'sem', 'figHandle', fig, ...
              'color_line', [r g b], 'color_area', [r g b])
```

`error` options: `'sem'`, `'std'`, `'c95'` (95% CI), `'c99'` (99% CI)

### `runPhotoBleachingRemoval` — drift correction

```matlab
signals = runPhotoBleachingRemoval(signals, 'lpCutOff', 0.1, ...
                                   'samplingRate', fs, 'filterOrder', 2)
```

Method options: `'lowpass'` (default), `'exp'` (exponential fit), `'spline'`

---

## Demo Files

All demos follow the same structure:

1. **Set paths** — edit `mainPath` / `folder` at the top of the script
2. **Load** — `loaduSMAART2mat` reads all `.mat` files in the folder and concatenates them
3. **Condition** — notch filter → bandpass → photobleaching removal
4. **Unmix** — `umxCONV` removes hemodynamic contamination
5. **Visualize** — PSD plots, raster/ERP plots, trigger-averaged traces

> `demo_1S1A` and `demo_2S2A` use identical code structure; the difference lies in which data columns are assigned to `sig` / `ref`. Adapt the column indices to match your hardware channel assignment.

---

## Known Limitations & Notes

- **Hard-coded sampling rates in demos**: the demos assume `fs` is read from `meta.fs`; verify this matches your acquisition settings.
- **FastICA dependency**: `umxICA` will error if FastICA is not on the MATLAB path. All other unmixing functions are self-contained.
- **`demo_1S1A` vs `demo_2S2A`**: these two demos currently share the same channel indexing template. Adapt `sig = signals(:, N)` and `ref = signals(:, M)` to reflect your actual configuration.
- **`umxICA` status**: this function is experimental. The unmixing coefficient derivation from the ICA mixing matrix may require tuning for specific datasets.
- **Intan reader**: `read_Intan_RHD2000_file.m` is a complete parser for the RHD2000 binary format; it is not modified from the original Intan distribution.

---

## License

GNU General Public License v3 — see [LICENSE](LICENSE).
