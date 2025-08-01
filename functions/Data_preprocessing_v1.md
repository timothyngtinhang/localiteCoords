# TMS Data Pre-processing Quality Control Manual

## Goal
To verify whether the targeted TMS coordinate matches the actual TMS stimulation site by organizing and analyzing relevant data files.

## Study Background
This manual covers the quality control procedures for a repeated-measures TMS study where:

- Each participant underwent TMS stimulation three times
- Each session targeted a different brain location
- The study measured impacts on cognitive task performance

## Challenges with TMS Data Files

### File Organization Issues
- The Localite TMS system generates numerous non-essential files during sessions
- Critical data files are scattered across various locations within each subject's directory
- Missing files are common due to technical or operator errors

### File Characteristics
- All relevant data files are in XML format
- Key file types include:
  - **EntryTarget**
  - **InstrumentMarker**
  - **TMS_Trigger**

### File Validity Criteria
- **EntryTarget** and **InstrumentMarker**: The last file created on a given date usually contains valid data
- **TMS_Trigger**: Files exceeding a specific size threshold typically contain complete data

## Summary of Goals and Challenges
1. **Quality Control**: Ensure the TMS coordinate targeted matches the actual stimulation site.
2. **Session Structure**: Each participant has three TMS sessions targeting different brain locations.
3. **Data Organization**: Large numbers of irrelevant files are generated, often stored in scattered subfolders.


4. **File Validation**: Valid files are generally the final entries (in the case of EntryTarget/InstrumentMarker) or those exceeding a size threshold (TMS_Trigger).

## Data Inspection

During data collection, target folder files may be saved in inconsistent or deeply nested directories. As an illustration:

```
data_localite/
└── 153/
    ├── 153_20240327_153_3bd8ee71f2d9093b/
    │   └── Sessions/
    │       └── Session_20240327172254667/
    │           └── InstrumentMarkers/
    │               └── InstrumentMarker20240327172502000.xml
    └── session1_20240202_153_68a4e7d0/
        ├── Sessions/
        │   └── Session_20240202134236974/
        │       ├── EntryTarget/
        │       │   └── EntryTarget20240202134412000.xml
        │       └── InstrumentMarkers/
        │           └── InstrumentMarker20240202134533000.xml
        └── session1_20240202_153_68a4e7d0/
            └── Sessions/
                └── Session_20240202164438019/
                    └── InstrumentMarkers/
                        └── InstrumentMarker20240202164555000.xml
```

These inconsistencies not only increase manual workload but also elevate the risk of file loss or misclassification. Variations in session naming and redundant subfolders further complicate inspection and analysis.

To address this, we can apply a folder-cleaning script that extracts only the relevant files and places them into a streamlined directory structure. The result is a clean, minimal, and navigable folder format:


```
Organized_Data/
└── 153/
    └── Sessions/
        └── Manually_created_sessions/
            ├── EntryTarget/
            │   ├── EntryTarget20240202134412000.xml
            │   └── EntryTarget20240311171732000.xml
            ├── InstrumentMarkers/
            │   ├── InstrumentMarker20240202134533000.xml
            │   ├── InstrumentMarker20240311171901000.xml
            │   ├── InstrumentMarker20240327172502000.xml
            │   └── InstrumentMarker20240202164555000.xml
            └── TMS_Trigger/
                └── (Only valid large-size files, if any)
```

### Folder-cleaning script
1. **Goal**: Create a coherent folder structure with only the three target subfolders: 
   - `EntryTarget`
   - `InstrumentMarker`
   - `TMS_Trigger`
2. **Method**: A script copies the original folder content while filtering out irrelevant files, placing only the valid target files into the respective subfolders.
3. **Result**: Easier inspection of missing files by date.

**Script Reference**: `first_extract_directory_tree_for_specified_folder.m`

