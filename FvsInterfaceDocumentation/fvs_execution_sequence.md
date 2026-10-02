# FVS Execution Sequence and Stop Points

## Overview

The Forest Vegetation Simulator (FVS) follows a structured execution sequence that can be interrupted at specific "stop points" to allow external programs to examine or modify data. This document describes the complete execution flow and identifies where each stop point occurs within the sequence.

## Main Execution Flow

### Initial Setup and Data Loading

FVS begins execution in the main `FVS()` subroutine, which serves as the primary entry point. The initialization phase consists of several critical steps:

**Command Line Processing**: FVS first processes command line parameters through `fvsSetCmdLine`, which parses options like `--keywordfile` and `--stoppoint`. If restart functionality is enabled, the system checks for existing restart files and sets appropriate restart codes.

**Restart Logic**: The `fvsRestart` function determines whether this is a fresh run or a continuation from a previous stop point. Restart codes can be positive (indicating the stop point to resume from) or negative (signaling that the calling program should make another FVS call to continue).

**Core Initialization**: `INITRE` initializes the FVS runtime environment, setting up internal data structures and reading configuration parameters.

**Tree Record Loading**: `NOTRE` loads and processes the initial tree inventory data, converting input records into FVS's internal tree representation. This is where the basic tree data (species, diameter, height, etc.) is read and validated.

### Stop Point 7: Post-Input, Pre-Imputation

**Stop Point 7** occurs immediately after `NOTRE` completes. At this point, FVS has loaded all tree records but has not yet performed any missing value imputation or model calibration. This stop point is specifically designed for external programs that need to examine the raw input data before FVS begins filling gaps or applying default values.

The description "after input read but before missing values are imputed" is precise - tree heights, crown ratios, and other variables that may be missing from the inventory have not yet been estimated or filled in.

### Model Calibration and Gap Filling

After resuming from Stop Point 7 (or continuing in normal execution), FVS performs several critical initialization steps:

**Root Disease Initialization**: `RDMN1` sets up the Western Root Disease Model if active.

**Dead Tree Processing**: `MPSDLP` and `DFBINV` compute dead trees per acre for various pest models.

**Calibration**: `CRATET` calibrates growth functions and fills gaps in the data. This is where missing tree heights, crown ratios, and other essential variables are imputed using FVS's built-in estimation procedures.

**Crown Width Calculation**: `CWIDTH` computes initial crown width values for all trees.

**Volume Calculations**: `VOLS` computes initial volume statistics for the stand.

### Statistical Processing and Initial Output

The system then computes various statistical summaries and generates initial output:

**Stand Statistics**: `STATS` computes statistical descriptions of the input data, and `DISPLY` writes the initial stand composition tables.

**Establishment Setup**: Various models are initialized for tree establishment, pest dynamics, and other ecological processes.

## Growth Cycle Processing

### Cycle Initiation

For each projection cycle (typically representing 5-10 years), FVS calls `TREGRO` (Tree Growth), which orchestrates the growth simulation process. `TREGRO` in turn calls `GRINCR` (Growth Increment), which handles the detailed calculation of individual tree growth and mortality.

### Growth Increment Calculation (GRINCR)

The `GRINCR` subroutine represents the heart of FVS's growth simulation and contains multiple stop points that allow detailed examination and modification of the growth process.

**Site Index Processing**: The routine first processes any site index modifications through `SETSITE` options, adjusting productivity parameters for the current cycle.

**Stand Density Calculations**: `SDICAL` and `SDICLS` compute Stand Density Index values, which are used throughout the growth calculations to assess competition levels.

### Stop Point 1: Pre-Event Monitor

**Stop Point 1** occurs just before the first call to the Event Monitor. At this point, all site parameters have been set, stand density has been calculated, but no management activities have been evaluated. This stop point allows external programs to examine or modify stand conditions before any management events are triggered.

The Event Monitor (`EVMON`) is FVS's system for triggering management activities based on stand conditions. Stop Point 1 provides the last opportunity to modify conditions before these triggers are evaluated.

### Stop Point 2: Post-Event Monitor Phase I

**Stop Point 2** occurs immediately after the first Event Monitor call. Any management activities triggered by the initial stand conditions have been scheduled, but not yet executed. This allows examination of what management activities will occur in this cycle.

### Management Activity Processing

Following the Event Monitor, FVS processes scheduled management activities through the `CUTS` subroutine, which handles thinning, harvesting, and other silvicultural treatments.

**Pre-Treatment Processing**: Economic evaluations (`ECSTATUS`) occur before treatments are applied.

**Treatment Application**: `CUTS` applies scheduled treatments, removing trees and modifying stand structure.

**Post-Treatment Updates**: Stand statistics are recalculated (`DENSE`) to reflect the new conditions after treatment.

### Stop Point 3: Post-Treatment

**Stop Point 3** occurs after all management activities have been applied but before the second phase of Event Monitor processing. This allows examination of stand conditions immediately following treatment application.

### Stop Point 4: Post-Event Monitor Phase II

**Stop Point 4** occurs after the second Event Monitor call, which can trigger additional management activities based on post-treatment conditions. This represents the final opportunity to examine conditions before growth calculations begin.

### Growth and Mortality Calculations

The core biological processes are then simulated through a series of specialized subroutines:

**Diameter Growth**: `DGDRIV` calculates diameter increment for each tree based on site conditions, competition, and species-specific growth patterns.

**Height Growth**: `HTGF` computes height increment, often as a function of diameter growth and site productivity.

**Regeneration**: `REGENT` handles the establishment of new trees, including both natural regeneration and planted seedlings.

**Mortality**: `MORTS` calculates the probability of death for each tree based on competition stress, age, and other factors.

### Stop Point 5: Post-Calculation, Pre-Application

**Stop Point 5** occurs after all growth and mortality calculations have been completed but before these changes are applied to the tree records. At this point, the growth increments (diameter and height) and mortality probabilities exist in temporary arrays but have not yet modified the actual tree data.

This stop point is particularly valuable for external growth models that want to examine or override FVS's calculated values before they are applied to the forest stand.

### Growth Application and Establishment

**Growth Application**: `GRADD` applies the calculated increments to tree diameters and heights, and removes trees that died during the cycle.

### Stop Point 6: Pre-Establishment

**Stop Point 6** occurs just before the establishment routines are called. All existing trees have been grown and mortality has been applied, but new trees have not yet been added to the stand. This allows examination of the stand after growth but before regeneration adds complexity.

**Tree Establishment**: Various establishment routines add new trees to the stand based on ecological conditions and management prescriptions.

## Cycle Completion and Output

After all biological processes are complete, FVS generates output for the current cycle:

**Statistical Updates**: Stand statistics are recalculated to reflect all changes that occurred during the cycle.

**Output Generation**: Various output tables and reports are generated, including stand composition tables, volume summaries, and specialized reports for active extensions.

**Visualization**: Stand visualization data is updated for graphical output systems.

## Projection Completion

This cycle repeats for each projection period until the specified end year is reached. Upon completion, FVS generates final summary reports and closes all output files.

## Stop Point Implementation Details

### File-Based vs. Memory-Based Stop Points

Stop points can be implemented in two ways:

**File-Based Stop Points**: Created using command line syntax like `--stoppoint=7,2025,restart.dat`. These create restart files that save the complete stand state, allowing FVS to be completely shut down and restarted later. File-based stop points generate negative restart codes (-1 through -7) that require two calls to FVS: the first returns the negative code, and the second performs the actual restart.

**Memory-Based Stop Points**: Created using the API function `fvsSetStoppointCodes()`. These maintain the stand state in memory and generate positive restart codes (1 through 7) that allow immediate continuation without file I/O.

### Practical Applications

Each stop point serves specific purposes in forest modeling applications:

- **Stop Point 7**: Data validation and preprocessing before FVS applies defaults
- **Stop Point 1**: Stand condition assessment before management evaluation  
- **Stop Points 2-4**: Management activity review and modification
- **Stop Point 5**: Growth model substitution or calibration
- **Stop Point 6**: Regeneration model replacement or enhancement

This structured approach allows external programs to integrate seamlessly with FVS's simulation process, enabling custom models, data validation, and specialized output generation while leveraging FVS's robust forest dynamics framework.