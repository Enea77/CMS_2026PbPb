# CMS 2026 PbPb Jet Health & Particle Flow Analysis

This repository contains an optimized ROOT framework for processing, monitoring, and plotting Jet Health parameters and Jet Particle Flow (PF) energy fractions in Lead-Lead (PbPb) heavy-ion collisions for the CMS experiment. 

The framework handles multi-dimensional correlation mapping via `THnSparseF` grids to isolate detector inefficiencies, compare different run periods (e.g., assessing the tracking performance with the Forward Pixel detector in vs. out), evaluate data-vs-simulation performance, and extract physics profiles as a function of event centrality.

---

## Workflow Architecture & Data Flow

The analysis is split into a robust **histogram production stage** and two multi-canvas **slide plotting modules**:

┌────────────────────────────────────────────────────────┐│           Input CMS Trees / Data Files                 ││ (hiEvtAnalyzer/HiTree, skimanalysis/HltTree, etc.)     │└───────────────────────────┬────────────────────────────┘│▼ [Executable Stage]┌────────────────────────────────────────────────────────┐│            JetHealth_PbPb_lxplus.cpp                   ││   - Applies vertex, event filter, & jet selections     ││   - Builds 'hjetkin' & 'hjetpf' THnSparseF objects     │└───────────────────────────┬────────────────────────────┘│▼ [Persistency Layer]┌────────────────────────────────────────────────────────┐│                     output.root                        │└───────────────┬────────────────────────┬───────────────┘│                        │▼ [Plotting Stage A]     ▼ [Plotting Stage B]┌──────────────────────────────┐  ┌──────────────────────────────┐│PlotJetHealthEtaPhiRegion.cpp │  │    PlotJetFractions.cpp      ││                              │  │                              ││ Generates side-by-side 2D    │  │ Generates 3-panel slides     ││ detector occupancy grids for │  │ showing Centrality-dependent ││ running comparisons (colz).  │  │ PF distributions & Ratios.   │└──────────────────────────────┘  └──────────────────────────────┘
### 1. Data Processing Loop (`JetHealth_PbPb_lxplus.cpp`)
* **Purpose**: Loops over events inside user-defined input file lists, checks vertex ranges ($|v_z| < 15\text{ cm}$), verifies noise/event filter branches (`ppvF`, `pclustF`, `pphfF`), ensures the Minimum Bias trigger fires, and enforces an analysis-level jet kinematic threshold ($p_T > 10\text{ GeV/c}$).
* **Outputs**: An organized `output.root` file containing multi-dimensional diagnostic data arrays.

### 2. Geometric Occupancy Diagnostics (`PlotJetHealthEtaPhiRegion.cpp`)
* **Purpose**: Compares spatial performance between two dataset contexts (e.g., tracking runs or subtraction configurations). It projects kinematic grids into $\eta$-$\phi$ space to visualize physical acceptance drops or masked regions.

### 3. Energy Fraction Profiles (`PlotJetFractions.cpp`)
* **Purpose**: Extracts detailed physics monitoring of jet compositions across distinct collision centrality blocks, producing publication-ready 3-panel slides showing comparative data trends along with their mathematical ratios.

---

## Core Framework Modifications & Technical Enhancements

This repository has been upgraded from its baseline iteration to transition from un-correlated kinematic plots into cross-correlated Particle Flow mapping.

### 1. High-Dimensional Extension of `THnSparseF` Arrays
The key optimization is a multi-dimensional array reconstruction inside `header/JetHealthHistograms.h`. The configuration of the sparse histograms is structured as follows:

* **Kinematic Sparse Histogram (`hjetkin`)**: A 4D sparse array tracking coordinate variables:
  $$\text{Axes:}\quad [0]\ p_T \quad\longrightarrow\quad [1]\ \eta \quad\longrightarrow\quad [2]\ \phi \quad\longrightarrow\quad [3]\ \text{hiBin (Centrality)}$$
* **Particle Flow Fractions Sparse Histogram (`hjetpf`)**: Expanded into an integrated **6D sparse matrix** designed to map jet constituent energy profiles:
  $$\text{Axes:}\quad [0]\ \text{pfFrac} \quad\longrightarrow\quad [1]\ \text{pfType} \quad\longrightarrow\quad [2]\ \eta \quad\longrightarrow\quad [3]\ \phi \quad\longrightarrow\quad [4]\ p_T \quad\longrightarrow\quad [5]\ \text{hiBin}$$

### 2. Hardcoded PF Index Token Maps
Particle Flow components are indexed via an internal type enumerator mapped across six dimensional bounds within `FillPF()` loops:

| Enumerator Index (`pfType`) | Acronym Shortname | Descriptive Axis Title |
| :--- | :--- | :--- |
| **`0`** | `CHF` | Charged Hadron Fraction |
| **`1`** | `NHF` | Neutral Hadron Fraction |
| **`2`** | `CEF` | Charged EM Fraction |
| **`3`** | `NEF` | Neutral EM Fraction |
| **`4`** | `MUF` | Muon Fraction |

### 3. Memory-Safe Automated Projection Pipelines
To abstract complex slicing configurations out of the core plotting scripts, `header/Utilities.h` implements two memory-resilient wrapper functions (`ProjectTHn1D` and `ProjectTHn2D`). 
* **Dynamic Range Intercepts**: These methods consume vector structures of `SparseRange` boundaries, iterate through active dimensions via `.SetRangeUser()`, apply a subtle upper-bound error offset adjustment (`hi - 0.0001`) to prevent floating-point bin leakage, extract a clone of the matrix projection, and immediately clear the state of the active axes via `.SetRange(0, -1)`. This architecture completely eliminates axis bleed-through bugs during downstream evaluation loops.

---

## Detailed Macro Explanations

### PlotJetHealthEtaPhiRegion.cpp
Designed to monitor localized performance drops across specific regions of the tracker or calorimeters.
* **Projections**: Slices `hjetkin` down to a 2D $\eta$-$\phi$ map (Axis 1 vs Axis 2), applying programmatic cuts to Jet $p_T$ (Axis 0) and Centrality `hiBin` (Axis 3).
* **Plot Structure**: Projects a dual-pane canvas (`1600 x 800`) splitting two target runs side-by-side using a uniform 2D `colz` representation.
* **Canvas Styling & Metadata Rules**: 
  * Normalizes histograms by their total integral count to cleanly compare shape variations.
  * Intercepts maximum/minimum $Z$-axis bounds across both histograms to enforce a synchronized color scale.
  * Draws normalized metadata layouts displaying `#bf{CMS} #it{Internal}`, the explicit jet algorithmic configuration (`akCs4PF`), kinematic limits ($p_T > \text{Cut}$), and centrality bounds inside the top margins.

### PlotJetFractions.cpp
Extacts physics profiles characterizing internal jet composition trends.
* **Projections & Conversion Pipelines**: 
  1. Filters `hjetpf` by selecting a chosen `pfType` interval on Axis 1, applying a minimum momentum threshold on Axis 4, and choosing a specific centrality range on Axis 5.
  2. Extracts a 2D histogram mapping Particle Flow fraction vs. Pseudorapidity $\eta$ (Axis 0 vs Axis 2).
  3. Converts the 2D plane into an average profile distribution by invoking `.ProfileX()`.
  4. Projects this average configuration into a clean 1D tracking representation using `.ProjectionX()`.
* **Plot Structure**: Generates a wide 3-panel visualization canvas (`3600 x 1200`):
  * **Panel 1**: Absolute Particle Flow distribution profile for Dataset A across multiple color-coded centrality layers (`hiBins`).
  * **Panel 2**: Absolute Particle Flow distribution profile for Dataset B under identical constraints.
  * **Panel 3**: A direct point-by-point overlay showing the mathematical ratio of Dataset A / Dataset B.
* **Visualization Annotations**: Automatically overlays vertical dotted lines at specific pseudorapidity reference points ($\eta = -2.3, -1.6, 1.5, 2.1$) to visually isolate and cross-reference geometric layout boundaries of the CMS subdetectors directly inside the physics slides.

---

## Technical Execution Guide

### 1. Running the Production Pipeline
The entry macro accepts configuration inputs natively via the command line or an interactive ROOT session:
```bash
# Compilation and execution via shell argument execution:
# Usage: ./JetHLT <filelist.txt> <output.root> <isMC>
root -l -b -q 'executable/JetHealth_PbPb_lxplus.cpp+("filelist.txt", "output_Data.root", false)'
2. Generating Comparative LayoutsExecute the visualization scripts to process local inputs and save generated diagnostics directly to timestamped directories:Bash# Run tracker occupancy comparisons
root -l executable/PlotJetHealthEtaPhiRegion.cpp

# Run jet fraction profile ratio tracking
root -l executable/PlotJetFractions.cpp

---
