clear, close all

%%%%%%%%%%%% Specify your configuration.m file here %%%%%%%%%%%%%%%%%%%%%%%
config_script = [pwd '\..\functions\config\Tim_TMS_pipeline_config.m']; % Change this line to switch configs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make config_script global so other functions can access it
global GLOBAL_CONFIG_SCRIPT
GLOBAL_CONFIG_SCRIPT = config_script;

config = pipeline_config();

% Ensure the current directory is the project root
cd(fileparts(mfilename('fullpath'))); % Change to the directory where run_all.m is located
cd('..'); % Go up one level to the project roots

%% This script transforms and extract raw TMS data into coordinates and give a view on the TMS data quality, by returning the euclidean distances 
% to add, total duration time between triggers start / end
% script to retrieve instrument marker by name not last file

%% Defining variables
scriptDir    = fileparts( mfilename('fullpath') );   % e.g. .../projectRoot/scripts
functionDir = fullfile(scriptDir, '..', 'functions');
configDir = fullfile(functionDir, 'config');

addpath(genpath( functionDir ) );  % go up, then into 'functions'
addpath(genpath( configDir ) );  % add config directory to path

% Create organized directory if it doesn't exist
if ~exist(fullfile(scriptDir, '..', config.dirs.organized), 'dir')
    mkdir(fullfile(scriptDir, '..', config.dirs.organized));
end

% Load TMS conditions and specifications from CSVs
tms_contrast = tms_conditions_config(config); 
condMap = config.input.tms_conditions_table;

rawDir = fullfile(scriptDir, '..', config.dirs.raw);  
organizedDir = fullfile(scriptDir, '..', config.dirs.organized); 

% rawDir = '../data/raw';
subjects = config.subjects;

%% you may choose what to run on each steps

% step 1: Create new organized Directory for all EntryTarget,
% InstrumentMarkers and TMSTrigger files

% organizeData(rawDir, organizedDir, subjects);

% step 2: return the filepaths of the organized directory at step 1
size_threshold = config.size_threshold; % TMSTrigger exceeding 100KB are selected
pathTable = extractPaths(organizedDir, subjects, size_threshold);

% ToDo: Instrumentmarker by names

% step 3: 
% Use the should-be template (when the session actually take place) to determine 
% the valid relevant file from the list of paths generated in step 2.
% in the following, it's from Google calendar

% finalizing the template here and should-be status here "scheduled",
% "Finished", "Not Scheduled", "Withdrawn".

% u need to input the valid_session_gcal.csv separately

SessPerSub = 3;
if config.utils.run_validSessWide2Long
    validTable_wide = readtable(fullfile(configDir, config.utils.validSess.wide_file), Delimiter=","); % u may skip this step if you have long data directly
    validTable_long = validSessWide2Long(validTable_wide, SessPerSub);
else
    validTable_long = readtable(fullfile(configDir, config.utils.validSess.long_file), Delimiter=","); 
end

% alternatively: validTable_long = readtable('valid_session_gcal_long.csv', Delimiter=",");

% condMap = table( ...
%     ["Vertex"; "iTBS"; "cTBS"], ...                            % Column #1: Condition
%     { ["EntryTarget","InstrumentMarkers","TMSTrigger"];        % Column #2: Types
%       ["InstrumentMarkers","TMSTrigger"];
%       ["InstrumentMarkers","TMSTrigger"] }, ...
%     'VariableNames', {'Condition','Types'} ...
% );

all_valid_files = mapSessionFiles(validTable_long, pathTable, condMap);

if config.utils.run_manualReplaceFile
    overwrite_csv = readtable(fullfile(configDir, config.utils.manualReplace.replacement_map_file), delimiter=",");
    all_valid_files = manualReplaceFile(all_valid_files, overwrite_csv);
end

all_files_with_coord = readCoordFromFiles(all_valid_files);
if config.utils.run_checkDurationMatch
   checkDurationMatch(all_files_with_coord, config);
end

% step 5: 

% optional, since iTBS and cTBS should share the same InstrumentMarker, but
% sometimes only one of them appears, insmarker_equivalence establish such
% equivalence

all_files_with_coord_eq = all_files_with_coord;
if config.utils.run_insmarker_equivalence
    all_files_with_coord_eq = insmarker_equivalence(all_files_with_coord, config.utils.insmarker_equivalence.markers);
end

% Run CsvCompareCoords function
if config.utils.run_CsvCompareCoords
    [CsvCompareCoords, CsvSummaryTable] = CsvCompareCoords(all_files_with_coord_eq, 'native_coords_summary.csv');
end

% every distance condition --------------------------
% tms_contrast = struct( ...
%     'GreenEntry', struct('RefType','EntryTarget',  'RowCond',"Vertex", 'idx',1 ), ... % Green Entry Point for Vertex condition
%     'RedTarget',  struct('RefType','EntryTarget',  'RowCond","Vertex", 'idx',2 ), ... % Red Target Point for Vertex condition
%     'InsEntry',   struct('RefType','InstrumentMarkers','RowCond',"Vertex", 'idx',1 ), ... % Instrument Marker Entry for Vertex
%     'OptiTBS',    struct('RefType','InstrumentMarkers','RowCond',"iTBS", 'idx',1 ), ...    % Optimized iTBS Instrument Marker
%     'OptcTBS',    struct('RefType','InstrumentMarkers','RowCond',"cTBS", 'idx',1 ));     % Optimized cTBS Instrument Marker

[TMS_Quality_raw, TMS_Quality_sum] = calculateCoords(all_files_with_coord_eq, tms_contrast);

% % Save extracted_paths.mat
save(fullfile(config.dirs.output_save_path, config.output.extracted_paths), 'pathTable');
% 
% % Save all_valid_files.mat
save(fullfile(config.dirs.output_save_path, config.output.valid_files), 'all_valid_files');
% 
% % Save all_files_with_coord.mat
save(fullfile(config.dirs.output_save_path, config.output.all_files), 'all_files_with_coord_eq');
% 
% % Save the final TMS_Quality table to a MAT file
save(fullfile(config.dirs.output_save_path, config.output.tms_quality_raw), 'TMS_Quality_raw');
save(fullfile(config.dirs.output_save_path, config.output.tms_quality_sum), 'TMS_Quality_sum');

fprintf('finished!\n');


function config = pipeline_config()
    global GLOBAL_CONFIG_SCRIPT
    
    if isempty(GLOBAL_CONFIG_SCRIPT)
        error('Global config script path not set. Run the main script first.');
    end
    
    run(GLOBAL_CONFIG_SCRIPT);
end
