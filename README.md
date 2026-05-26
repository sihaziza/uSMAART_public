# uSMAART_public
To process and analyze fiber photometry data using uSMAART (ultra‑Sensitive Measurement of Aggregate Activity in Restricted cell‑Types). uSMAART is a laser‑based fiber‑photometry optical instrument specifically designed to detect high‑frequency voltage dynamics from identified, user-targeted cell-types using genetically-encoded fluorescent voltage indicators. uSMAART can measure up to 2 fluorescent sensor signals from up to two brain areas at once. 
This repo provides MATLAB code to preprocess, analyze, and visualize uSMAART recordings across common acquisition configurations.

Demo files highlighting processing steps for each situation are provided (1S1A, 2S1A, and 2S2A, for S=sensor and A=brain area).

# Quick start

git clone https://github.com/sihaziza/uSMAART_public.git
cd uSMAART_public
% add the repo to your MATLAB path (e.g., savepath after "addpath(genpath(pwd))")

Run a demo in MATLAB:

% Single sensor, single area (1S1A)
demo_1S1A

% Two sensors, one area (2S1A)
demo_2S1A

% Two sensors, two areas (2S2A)
demo_2S2A

% Imaging demo
demo_imaging

# Citing uSMAART

If this code contributes to your work, please cite:

Haziza S., et al. Imaging high-frequency voltage dynamics in multiple neuron classes of behaving mammals. Cell, 2025. (https://doi.org/10.1016/j.cell.2025.06.028)
