function pathTable = extractPaths(organizedDir, subjectIDs, sizeThreshold)
% EXTRACTPATHS Extract file paths from organized TMS data
%
% Syntax:
%   pathTable = extractPaths(organizedDir, subjectIDs, sizeThreshold)
%
% Inputs:
%   organizedDir  - String, path to organized data directory
%   subjectIDs    - Array, subject IDs to process
%   sizeThreshold - Numeric, minimum file size for TMSTrigger files
%
% Outputs:
%   pathTable     - Table with columns: Subject, Type, Date, Filename, FullPath
%
% Example:
%   paths = extractPaths('data/organized/', 401:441, 100000000);

%% fetch subj id
% List folders in the data directory
d = dir(organizedDir);

% Target folder types to search within each subject folder
targetFolderTypes = {'EntryTarget', 'InstrumentMarkers', 'TMSTrigger'};

% Minimum file size threshold for TMSTrigger files (in bytes)
MIN_FILE_SIZE_THRESHOLD = sizeThreshold;

%% INITIALIZE RESULTS STORAGE
% Arrays to store results
results = struct();
results.dates = {};
results.filenames = {};
results.paths = {};
results.subjects = {};
results.folderTypes = {};

%% MAIN PROCESSING LOOP
disp('Starting data processing...');
disp('=========================');

% Process each subject
for subjectIdx = 1:length(subjectIDs)
    currentSubjectID = subjectIDs(subjectIdx);
    subjectCode = ['CC', num2str(currentSubjectID)];
    
    disp(['Processing subject: ', subjectCode]);
    disp('-------------------------');
    
    % Process each folder type for this subject
    for folderTypeIdx = 1:length(targetFolderTypes)
        currentFolderType = targetFolderTypes{folderTypeIdx};
        disp(['  Processing folder type: ', currentFolderType]);
        
        % Find all folders for this subject, return a list of those folders
        subjectFolder = fullfile(organizedDir, num2str(currentSubjectID));
        subFolderPath = fullfile(subjectFolder, 'Sessions', 'Manually_created_sessions', currentFolderType);

        xmlFiles = dir(fullfile(subFolderPath, '*.xml'));

        names = {xmlFiles.name};

        
        if strcmp(currentFolderType, 'EntryTarget')

        elseif strcmp(currentFolderType, 'InstrumentMarkers')

        elseif strcmp(currentFolderType, 'TMSTrigger')
            filtered_name_bool = startsWith(names, 'TriggerMarkers_Coil0');
            names = names(filtered_name_bool);
        end

        datesFromNames = cellfun(@(x) x(end-20:end-13), names, 'UniformOutput', false); % the last few characters are always datestr

        uniqueDates = unique(datesFromNames);

        %% Case 1: find the large file in TriggerMarker_Coil0
        if strcmp(currentFolderType, 'TMSTrigger')
            for i=1:length(uniqueDates)
                currentUniqueDate = uniqueDates(i);
                ContainCurrentDate_bool = contains(names, currentUniqueDate);
                currentUniqueDateFiles = names(ContainCurrentDate_bool); % files matching the current unique date

                for j = 1:length(currentUniqueDateFiles) %% get the file info for all data in a single date
                    fullFilePath = fullfile(subFolderPath, currentUniqueDateFiles{j});
                    fileInfo = dir(fullFilePath);
                    
                    if fileInfo.bytes > 100000 %% filter out files with small sizes
                        % This file is larger than 100KB
                        fprintf('Large file found: %s (%d bytes)\n', currentUniqueDateFiles{j}, fileInfo.bytes);
                        
                        results.dates{end+1} = currentUniqueDate{1};
                        results.filenames{end+1} = currentUniqueDateFiles{j};
                        results.paths{end+1} = fullFilePath;
                        results.subjects{end+1} = subjectCode;
                        results.folderTypes{end+1} = currentFolderType;

                    end
                end
            end
        else %% Case 2: find the last created file in EntryTarget & InstrumentMarker
            for i=1:length(uniqueDates)
                currentUniqueDate = uniqueDates(i);
                ContainCurrentDate_bool = contains(names, currentUniqueDate);
                currentUniqueDateFiles = names(ContainCurrentDate_bool); % files matching the current unique date
                % initialize timeArray
                timeArray = zeros(length(currentUniqueDateFiles), 1);

                for j = 1:length(currentUniqueDateFiles)
                    fullFilePath = fullfile(subFolderPath, currentUniqueDateFiles{j});
                    fileInfo = dir(fullFilePath);

                    if ischar(fileInfo.name) 
                        TimeFromNames_cell = {fileInfo.name(end-12:end-4)};
                    else
                        TimeFromNames_cell = cellfun(@(x) x(end-12:end-4), fileInfo.name, 'UniformOutput', false); % the last few characters are always datestr
                    end
                    timeArray(j) = str2num(cell2mat(TimeFromNames_cell));    % convert to num by to mat and to int
                end
                [maxTimeValue, maxTimeIndex] = max(timeArray);
                latestFile = currentUniqueDateFiles{maxTimeIndex};
                fprintf('Latest file found on %s: %s\n', currentUniqueDate{1}, latestFile);

                results.dates{end+1} = currentUniqueDate{1};
                results.filenames{end+1} = latestFile;

                LatestFilePath = fullfile(subFolderPath, latestFile);

                results.paths{end+1} = LatestFilePath;
                results.subjects{end+1} = subjectCode;
                results.folderTypes{end+1} = currentFolderType;
            end
        end
    end
end

%% CREATE OUTPUT TABLE
pathTable = table(...
    results.subjects', ...
    results.folderTypes', ...
    results.dates', ...
    results.filenames', ...
    results.paths', ...
    'VariableNames', {'Subject', 'Type', 'Date', 'Filename', 'FullPath'});

% Save to CSV
config = pipeline_config();

outDir = fullfile('..','data','organized', config.dataset_name);
if ~exist(outDir,'dir')
    mkdir(outDir);
end

outputFileName = 'step_1_all_extracted_path.csv';
fullpath = fullfile(outDir, outputFileName);
writetable(pathTable, fullpath);
disp(['Results saved to ', outputFileName]);