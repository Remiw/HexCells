# HexCells

ROOT-based signal processing pipeline for the time-resolution analysis of hexagonal plastic scintillator detectors with Silicon Photomultiplier (SiPM) readout.

This repository contains the analysis framework developed for the beam-test studies presented in my undergraduate thesis, **"Time Resolution and Beam Test of Hexagonal Plastic Scintillator Counters with Silicon Photomultiplier Readout"**, carried out within the SND@HL-LHC collaboration.

The code converts raw oscilloscope waveforms into ROOT TTrees, reconstructs detector signals, extracts pulse parameters, determines arrival times using both Fixed Threshold and digital Constant Fraction Discriminator (dCFD) algorithms, and evaluates the detector time resolution from beam-test data.

The software was developed to characterize prototype hexagonal plastic scintillator counters for the proposed fast timing detector of SND@HL-LHC.

---

## Repository Structure

| File | Description |
|------|-------------|
| **rootfyTwoCh.C** | Converts raw waveform data into ROOT TTrees. |
| **WaveFormFunctions.C** | Waveform processing library including baseline estimation, pulse reconstruction, amplitude, charge, pulse width, Fixed Threshold and dCFD timing algorithms. |
| **mergedata.C** | Main analysis macro for detector characterization and time-resolution studies. |
| **Treed1.root, Treed2.root, Treed3.root** | Example ROOT files provided to reproduce the analysis. |
| **try.root** | Additional test ROOT file. |

---

## Analysis Workflow

```
Raw Oscilloscope Data
        │
        ▼
 rootfyTwoCh.C
        │
        ▼
    ROOT TTrees
        │
        ▼
WaveFormFunctions.C
        │
        ▼
 Pulse Reconstruction
 Baseline Correction
 Pulse Parameters
 Fixed Threshold
 dCFD Timing
        │
        ▼
 mergedata.C
        │
        ▼
 Time Resolution Analysis
```

---

## Software

- CERN ROOT Framework (C++)
- ROOT TTrees

---

## Note

This repository currently contains only the beam-test analysis developed for the thesis.

The cosmic-ray analysis is maintained in a separate repository and will be integrated into this project in a future update to provide the complete analysis framework.

---

## Context

This work was developed as part of the detector characterization activities supporting the fast timing detector proposed for **SND@HL-LHC**.
