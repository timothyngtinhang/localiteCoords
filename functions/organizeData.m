function organizeData(rawDir, organizedDir, validSubjects)
% organizeData  Copy XML files from rawDir into organizedDir by subject & type.
%
%   organizeData(rawDir, organizedDir, validSubjects) does the following:
%     • rawDir          – root folder containing subject subfolders
%     • organizedDir    – target base folder where “Organized_Data” will be created
%     • validSubjects   – numeric vector of allowed subject IDs (e.g. [401:441])
%
%   If you omit any inputs, defaults are:
%     rawDir       = pwd
%     organizedDir = fullfile(pwd, 'Organized_Data')
%     validSubjects = 401:410

    %-----------------------------
    % 1) Handle missing inputs
    %-----------------------------
    % if nargin < 1 || isempty(rawDir)
    % 
    % end
    % if nargin < 2 || isempty(organizedDir)
    %     organizedDir = fullfile(rawDir, 'Organized_Data');
    % end

    if nargin < 3 || isempty(validSubjects)
        validSubjects = 401:410;
    end
    %-----------------------------
    % 2) Prepare the target folder
    %-----------------------------
    if exist(organizedDir, 'dir')
        rmdir(organizedDir, 's');
    end
    mkdir(organizedDir);

    % %-----------------------------
    % % 3) Find subject folders in rawDir
    % %   (only those whose names are numeric and in validSubjects)
    % %-----------------------------
    allItems   = dir(rawDir);
    isFolder   = [allItems.isdir];
    names      = {allItems.name};
    % Exclude “.”, “..”, and the target folder (in case it lives under rawDir)
    excludeIdx = ismember(names, {'.','..',    ...
                                  fullfile('.', rawDir)}); 
    candidate  = allItems(isFolder & ~excludeIdx);

    % Convert each folder name to a number (NaN if not numeric)
    folderNums = cellfun(@(s) str2double(s), {candidate.name});
    validMask  = ismember(folderNums, validSubjects);
    subjectDirs = candidate(validMask);

    % Loop through each subject directory
    for i = 1:length(subjectDirs)
        subjName = subjectDirs(i).name;
        subjDir = fullfile(rawDir, subjName);
        fprintf('Processing subject folder: %s\n', subjName);
        
        % Create the new directory structure for this subject:
        % Organized_Data/subject_id/Sessions/Manually_created_sessions/...
        newSubjDir   = fullfile(organizedDir, subjName);
        sessionsDir  = fullfile(newSubjDir, 'Sessions');
        manualDir    = fullfile(sessionsDir, 'Manually_created_sessions');
        
        % Create new directories for each folder type
        newEntryTargetDir       = fullfile(manualDir, 'EntryTarget');
        newInstrumentMarkersDir = fullfile(manualDir, 'InstrumentMarkers');
        newTMSTriggerDir        = fullfile(manualDir, 'TMSTrigger');
        
        if ~exist(newEntryTargetDir, 'dir'),       mkdir(newEntryTargetDir);       end
        if ~exist(newInstrumentMarkersDir, 'dir'), mkdir(newInstrumentMarkersDir); end
        if ~exist(newTMSTriggerDir, 'dir'),        mkdir(newTMSTriggerDir);        end
        
        % Define the folder types to search for
        folderTypes = {'EntryTarget', 'InstrumentMarkers', 'TMSTrigger'};
        
        % Loop over each folder type
        for f = 1:length(folderTypes)
            folderType = folderTypes{f};
            fprintf(' Searching for %s folders...\n', folderType);
            
            % Recursively get all items matching the folder type name.
            folderDirs = dir(fullfile(subjDir, '**', folderType));
            % Filter to get only actual directories.
            folderDirs = folderDirs([folderDirs.isdir]);
            
            if isempty(folderDirs)
                fprintf('  No %s folder found for subject %s\n', folderType, subjName);
            else
                % Loop over each found folder of the current folder type.
                for j = 1:length(folderDirs)
                    currentFolder = fullfile(folderDirs(j).folder, folderDirs(j).name);
                    fprintf('  Found %s folder: %s\n', folderType, currentFolder);
                    
                    % List all XML files with names starting with the folder type.
                    xmlFiles = dir(fullfile(currentFolder, [folderType(1) '*.xml'])); % match with first character of the foldertype (e.g. T*.xml for TMSTrigger folder)
                    if isempty(xmlFiles)
                        fprintf('    No XML files found in %s\n', currentFolder);
                    else
                        % Determine the new target directory for the current folder type.
                        switch folderType
                            case 'EntryTarget'
                                newTargetDir = newEntryTargetDir;
                            case 'InstrumentMarkers'
                                newTargetDir = newInstrumentMarkersDir;
                            case 'TMSTrigger'
                                newTargetDir = newTMSTriggerDir;
                        end
                        
                        % Loop through each XML file and copy it.
                        for k = 1:length(xmlFiles)
                            sourceFile  = fullfile(xmlFiles(k).folder, xmlFiles(k).name);
                            newFilePath = fullfile(newTargetDir, xmlFiles(k).name);
                            
                            copyfile(sourceFile, newFilePath);
                            fprintf('    Copied %s to %s\n', sourceFile, newFilePath);
                        end
                    end
                end
            end
        end
    end
    fprintf('finish executing organizeData.m\n')
end