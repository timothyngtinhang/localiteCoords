function config_output = tms_conditions_config(config)
% TMS_CONDITIONS_CONFIG defines the experimental conditions and coordinate specifications.
%
% This file centralizes the configuration for the TMS data processing pipeline.
% New users can easily modify experimental conditions, marker types, and
% coordinate comparison specifications without altering the main processing scripts.
%
% User Guide:

% 2.  **Defining Coordinate Specifications (condSpec):**
%     - The coordinate specifications are now loaded from `tms_condspec.csv`.
%     - Open `tms_condspec.csv` in a spreadsheet editor (e.g., Excel) to modify.
%     - Column 'ComparisonName': A descriptive name for the comparison (e.g., "GreenEntry").
%     - Column 'ConditionName': The experimental condition this applies to (must match a 'Condition' in `tms_conditions.csv`).
%     - Column 'ReferenceCoord': The type of marker to reference (e.g., 'EntryTarget', 'InstrumentMarkers').
%     - Column 'TargetCoord': The type of marker to target (e.g., 'TMSTrigger').
%     - Column 'EntryTarget': For 'EntryTarget' `ReferenceCoord`, specify 'Entry' or 'Target'. Use 'NA' otherwise.
%     - Example row in CSV:
%       `ComparisonName,ConditionName,ReferenceCoord,TargetCoord,EntryTarget`
%       `GreenEntry,Vertex,EntryTarget,TMSTrigger,Entry`
%       `OptiTBS,iTBS,InstrumentMarkers,TMSTrigger,NA`
%
% --------------------------------------------------------------------------------------------------

% read table
    
    csv_spec = readtable(fullfile(config.dirs.functions, 'config', config.input.tms_condspec_file), 'Delimiter', ',');
    
    % loop over each row and assign

    condSpec = struct();
    for i = 1:height(csv_spec)
        contrast = csv_spec.ComparisonName{i};
        cond = csv_spec.ConditionName{i};
        ref = csv_spec.ReferenceCoord{i};
        tar = csv_spec.TargetCoord{i};
        entrytarget = csv_spec.EntryTarget{i};

        idx = 1 + (strcmp(entrytarget, 'Target') && strcmp(ref, 'EntryTarget'));

        if ~strcmp(tar, 'TMSTrigger')
            fprintf('TargetCoord not specified as TMSTrigger. Will use TMSTrigger as target Coordinate nonetheless');
        end

        condSpec.(contrast) = struct( ...
            'RefType', ref, ...
            'RowCond', cond, ...
            'idx', idx);
    end

    config_output = condSpec;
end
