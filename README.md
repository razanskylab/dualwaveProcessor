# Processing and Analysis Code for Blood sO2 Datasets Acquired with Ultrasound-aided Spectroscopic Large-scale Optoacoustic Microscopy

## Initialization
Go to the cloned repo directory, and run startup to intialize the environment.

## 2D and 3D unmixing
- Run test_unmix2D to perform 2D spectral unmixing, i.e. just take maximum signal amplitude without considering depth- and wavelength-dependent fluence variations.
- Run test_unmix3D to perform 3D spectral unmixing, i.e. performs 1D fluence modeling and correction based on skin surface and vessel depth information, and unmix based on fluence-corrected maximum signal amplitudes.

## Vessel parameter quantification
After obtaining good sO2 maps, follow the sequence of scripts (in paperScripts folder) below to perform step-by-step vascular parameter quantification:
1. crop_wrap_sig, crops out plastic kitchen wrap signal in US dataset, as a preparation step for skin surface estimation;
2. estim_skin_surface, estimates the skin surface based on US dataset by identifying the first prominent reflection layer on a per B-scan basis;
3. segment_layered_vessels, segments vessels from different skin layers (e.g. epidermis, dermis, hypodermis) by overlaying the estimated skin surface onto the OA dataset;
4. get_wound_boundary, delineates the wound boundary in Fiji (ImageJ) by manually drawing a boundary of the wound center, indicated by a lack of vessel signals. Also generates concentric bands around the wound boundary for band-wise analysis.
5. get_structural_maps, extracts structural parameters from vessel maps using the open-source PostProGUI: https://github.com/razanskylab/PostProGUI;
6. quantify_wound_area_reduction, quantifies the reduction in wound areas over time by counting pixels inside wound boundaries;
7. illustrate_healing_front_granulation_area, generates figures to illustrate the definitions of healing front and granulation tissue area;
8. get_bandwise_stat, extracts statistics on so2 and structural vessel parameters per-band;
9. analyze_granulation_tissue, extracts statistics on sO2 and structural vessel parameters in the granulation tissue area over time;
10. generateBoxChart, generates box plots presented in the paper, based on the mat files containing extracted vessel parameters in matData folder.

NOTE: the script in step 10 can be run directly. The scripts for steps 1 to 9 and the two test scripts above require raw datasets that are too large to share here, but will be available upon request: weiye.li@uzh.ch.
