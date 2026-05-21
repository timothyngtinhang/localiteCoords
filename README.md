# Localite TMS Post-processing Workflow

## Description

This repository contains a MATLAB workflow for post-processing Localite TMS (Transcranial Magnetic Stimulation) data. It was developed to support research data checking by extracting coordinate information from Localite export files, mapping files to experimental sessions and conditions, calculating customized Euclidean-distance comparisons, and checking whether pulse durations match the expected values for each condition.

The workflow handles multiple Localite-related file types, including EntryTarget, InstrumentMarkers, and TMSTrigger files. It is mainly intended as a reference implementation for the original research workflow rather than a general-purpose, actively maintained package.

## Status

**Legacy project / no longer actively maintained.**

The core workflow was developed for a specific TMS research project and folder structure. It may require adaptation before reuse with other Localite exports, experimental designs, or data organization formats.

This repository is kept public as a reference for the original workflow and as an example of research data post-processing automation. Bugs may exist, because apparently research data pipelines are legally required to contain at least one cursed edge case.

## Features

- Processes multiple Localite/TMS file types:
  - EntryTarget
  - InstrumentMarkers
  - TMSTrigger
- Extracts 3D coordinate information from Localite-related files
- Maps files to subjects, sessions, and experimental conditions
- Supports customized coordinate comparisons across conditions
- Calculates Euclidean distances between reference coordinates and TMS trigger locations
- Checks whether pulse duration matches expected condition values
  - iTBS: approximately 188 seconds
  - cTBS: approximately 40 seconds
- Supports multiple experimental conditions, including Vertex, iTBS, and cTBS
- Uses CSV-based session validation and condition specification files
- Provides two execution modes:
  - Batch processing through `run_all.m`
  - Self-contained processing through `one_script.m`

## Installation

Clone the repository:

```bash
git clone https://github.com/timothyngtinhang/localiteCoords.git
cd localiteCoords
```

Ensure MATLAB is installed.

No additional MATLAB toolboxes are required.

**Quick Start**
---------------

### **Option 1:**

### **`run_all.m`**

This is the recommended option for processing multiple datasets or using an external configuration file.

1.  Edit the pipeline configuration file:
2.  Configure the main settings:
3.  Prepare the session validation file in `functions/config/`.

    Supported formats include:

    -   `Tim_optional_valid_session_list_wide.csv`
    -   `Tim_valid_session_list_long.csv`
4.  Run the pipeline:

```
cd script
run_all
```

### **Option 2:**

### **`one_script.m`**

This option keeps the main parameters and processing logic in one script. It is easier to inspect, but less modular and less reproducible than the configuration-based workflow.

1.  Edit the parameters directly inside `one_script.m`, including:
    -   Raw data directory
    -   Organized output directory
    -   Subject IDs
    -   Session validation file path
    -   Condition mapping
2.  Run the script:

```
cd script
one_script
```

**Expected Data Structure**
---------------------------

```
data/
├── raw/
│   └── [dataset_name]/
│       └── [subject_folders]/
└── organized/
    └── [dataset_name]/
        ├── step_1_all_extracted_path.mat
        ├── step_2_all_valid_files.mat
        ├── step_3_all_files_with_coord.mat
        ├── TMS_Quality_Raw.mat
        └── TMS_Quality_Summary.mat
```

**Processing Workflow**
-----------------------

```
flowchart TD
    A[Raw Localite Data<br/>EntryTarget, InstrumentMarkers, TMSTrigger] --> B{Choose Processing Mode}

    B -->|Batch Processing| C[run_all.m]
    B -->|Single Dataset| D[one_script.m]

    C --> E[Load Config<br/>Tim_TMS_pipeline_config.m]
    D --> F[Inline Parameters<br/>subjects, paths, conditions]

    E --> G[Step 1: organizeData<br/>Restructure raw files]
    F --> G

    G --> H[Step 2: extractPaths<br/>Find valid file paths<br/>Filter by size threshold]

    H --> I[Step 3: Session Validation<br/>validSessWide2Long<br/>Match with session schedule]

    I --> J[Step 4: mapSessionFiles<br/>Map files to conditions<br/>Vertex, iTBS, cTBS]

    J --> K[Optional: manualReplaceFile<br/>Override specific files]

    K --> L[Step 5: readCoordFromFiles<br/>Extract 3D coordinates]

    L --> M[Optional: insmarker_equivalence<br/>Share markers between conditions]

    M --> N[Step 6: calculateCoords<br/>Compute Euclidean distances<br/>Reference to target points]

    N --> P[Step 7: checkDurationMatch<br/>Verify pulse duration matches expected condition time]

    P --> O[Output Results<br/>TMS_Quality_Raw.mat<br/>TMS_Quality_Summary.mat]

    style C fill:#e1f5fe
    style D fill:#fff3e0
    style O fill:#e8f5e8
```

**Key Functions**
-----------------

### **Core Processing Functions**

-   `organizeData.m`\
    Restructures raw Localite files into an organized directory structure.
-   `extractPaths.m`\
    Scans the organized directory and returns file paths that meet file size and naming criteria.
-   `mapSessionFiles.m`\
    Maps extracted files to subjects, sessions, and experimental conditions.
-   `readCoordFromFiles.m`\
    Extracts 3D coordinate information from Localite-related files.
-   `calculateCoords.m`\
    Computes Euclidean distances between configured reference and target coordinates.
-   `checkDurationMatch.m`\
    Checks whether pulse durations match expected values for each TMS condition.

### **Utility Functions**

-   `validSessWide2Long.m`\
    Converts a wide-format session validation table into long format.
-   `manualReplaceFile.m`\
    Allows specific files to be manually replaced or overridden in the pipeline.
-   `insmarker_equivalence.m`\
    Defines equivalence between instrument markers across conditions when needed.

### **Configuration Files**

-   `Tim_TMS_pipeline_config.m`\
    Main configuration file for subjects, paths, file thresholds, and processing options.
-   `tms_conditions_config.m`\
    Loads experimental condition settings and coordinate specifications from CSV files.

**Configuration Example**
-------------------------

Example settings from `functions/config/Tim_TMS_pipeline_config.m`:

```
config.subjects = 401:441;
config.size_threshold = 100000;
config.dirs.raw = 'data/raw/Tim_TMS';
config.dirs.organized = 'data/organized/Tim_TMS';
config.run_organizeData = false;
```

**Session Validation Files**
----------------------------

The workflow supports session validation files to control which subject-session combinations should be included.

Examples:

-   `Tim_optional_valid_session_list_wide.csv`\
    One row per subject, with sessions represented across columns.
-   `Tim_valid_session_list_long.csv`\
    One row per subject-session combination.

**Condition Specification**
---------------------------

Condition-specific coordinate comparisons can be configured through CSV files such as:

-   `Tim_condspec.csv`

This file defines which coordinate sources should be compared for each condition, such as comparing EntryTarget coordinates with TMSTrigger coordinates.

**Output Files**
----------------

The main output files are:

-   `TMS_Quality_Raw.mat`\
    Detailed distance calculations for individual triggers or files.
-   `TMS_Quality_Summary.mat`\
    Aggregated quality metrics by subject, session, or condition.
-   `step_1_*.mat`, `step_2_*.mat`, `step_3_*.mat`\
    Intermediate files for checking and debugging the pipeline.

**Experimental Conditions**
---------------------------

The original workflow was designed around the following TMS-related conditions:

-   **Vertex**\
    Control condition using EntryTarget, InstrumentMarkers, and TMSTrigger files.
-   **iTBS**\
    Intermittent theta burst stimulation condition using InstrumentMarkers and TMSTrigger files.\
    Expected pulse duration: approximately 188 seconds.
-   **cTBS**\
    Continuous theta burst stimulation condition using InstrumentMarkers and TMSTrigger files.\
    Expected pulse duration: approximately 40 seconds.

**Troubleshooting**
-------------------

Common issues include:

1.  **TMSTrigger files not detected**\
    Adjust `config.size_threshold` if relevant files are being filtered out.
2.  **Missing coordinates**\
    Check whether the original Localite files contain valid coordinate data.
3.  **Session mismatch**\
    Ensure the session validation CSV matches the actual experimental timeline and subject IDs.
4.  **Path errors**\
    Use absolute paths or confirm that MATLAB is running from the expected working directory.
5.  **Pulse duration mismatch**\
    Check the output from `checkDurationMatch.m`, then manually inspect the corresponding raw TMSTrigger files if needed.

**Notes on Reuse**
------------------

This workflow was written for a specific research context. Before applying it to a new dataset, users should review and adapt:

-   Folder structure assumptions
-   Subject ID format
-   Session validation files
-   Condition labels
-   Coordinate comparison rules
-   Expected pulse durations
-   Localite export format

The code is provided as a reference workflow rather than a plug-and-play package.

