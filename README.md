# Kepler-1976 b Transit Analysis Pipeline

This repository hosts a comprehensive Python pipeline developed to process and analyze light curves from the **Kepler Space Telescope**. The project successfully extracts the transit signal of the exoplanet **Kepler-1976 b**, performs precise phase-folding, and estimates its planetary radius.

## Scientific Objectives
The project was designed with the following goals in mind:
* **Signal Isolation:** Removing long-term stellar variability and instrumental noise using 3rd-degree polynomial detrending.
* **Phase Synchronization:** Implementing an automated centering logic to align the transit event at orbital phase zero.
* **Radius Characterization:** Calculating the absolute radius of the planet in Earth units ($R_\oplus$) to determine its classification.

## Key Results & Analysis
The pipeline successfully identified a transit signal consistent with a **Gas Giant**:
* **Calculated Planet Radius:** $11.69~R_\oplus$ (similar to Jupiter).
* **Stellar Radius:** $118.05~R_\oplus$.
* **Light Curve Morphology:** The observed "V-shape" in the binned data highlights physical effects such as **Limb Darkening** and Kepler's 30-minute integration smearing.

## Features
* **Automated Centering:** Finds the minimum flux phase and applies a wrap-around shift to prevent "split" transits.
* **Scalability:** The pipeline is a generalized framework. By updating $P$, $T_0$, and $R_*$, it can characterize any target in the Kepler archive.
* **Robust Statistics:** Uses median binning to suppress photon noise and ensure robustness against outliers.

## Project Structure
* `kepler_analysis_narges.py`: The core Python script containing the processing pipeline.
* `Kepler1976b_Transit_Analysis.pdf`: Detailed scientific report of the findings.
* `ff5.png`: Visualization of the final results.
