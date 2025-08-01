
% PIPELINE_CONFIG Central configuration for TMS preprocessing pipeline

% Dataset settings
config.dataset_name = 'Tim_TMS'; % Name of the dataset to process

config.subjects = 401:441;

% search rule for TMSTrigger
config.size_threshold = 100000000; % 100KB in bytes
% search rule for InstrumentMarker and EntryTarget
% // find the last file on that specific date

config.coordinate_systems = {'RAS', 'LAS'}; % doesn't support for MNI
config.file_types = {'EntryTarget', 'InstrumentMarkers', 'TMSTrigger'};


% Directory paths (now dynamically set based on dataset_name)
config.dirs.raw = fullfile('data', 'raw', config.dataset_name);
config.dirs.organized = fullfile('data', 'organized', config.dataset_name);
config.dirs.functions = 'functions/';

% input file names (as contained in config directory)
config.input.tms_conditions_table = table( ...
                                    ["Vertex"; "iTBS"; "cTBS"], ...                            % Column #1: Condition
                                    { ["EntryTarget","InstrumentMarkers","TMSTrigger"];        % Column #2: File Types included in each condition
                                      ["InstrumentMarkers","TMSTrigger"];
                                      ["InstrumentMarkers","TMSTrigger"] }, ...
                                    'VariableNames', {'Condition','Types'} ...
                                );

config.input.tms_condspec_file = 'Tim_condspec.csv';

%% Utility script settings (optional)
config.run_organizeData = false; % Create simplified Directory structure (necessary for first run)
config.utils.run_validSessWide2Long = true; % Set to true to run validSessWide2Long.m
config.utils.run_manualReplaceFile = false;  % Set to true to run manualReplaceFile.m
config.utils.run_insmarker_equivalence = true; % Set to true to run insmarker_equivalence.m

% Parameters for validSessWide2Long.m
config.utils.validSess.wide_file = 'Tim_optional_valid_session_list_wide.csv';
config.utils.validSess.long_file = 'Tim_valid_session_list_long.csv';

% Parameters for manualReplaceFile.m
config.utils.manualReplace.replacement_map_file = 'Tim_optional_overwriting_valid.csv';

% Parameters for insmarker_equivalence.m
config.utils.insmarker_equivalence.markers = ["iTBS", "cTBS"];

%% Output file names
config.output.extracted_paths = 'step_1_all_extracted_path.mat';
config.output.valid_files = 'step_2_all_valid_files.mat';
config.output.all_files = 'step_3_all_files_with_coord.mat';
config.output.tms_quality_raw = 'TMS_Quality_Raw.mat'; 
config.output.tms_quality_sum = 'TMS_Quality_Summary.mat'; % New output for quality results

