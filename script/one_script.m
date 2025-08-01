clear, close all
%% This script transforms and extract raw TMS data into coordinates and give a view on the TMS data quality, by returning the euclidean distances 
% to add, total duration time between triggers start / end
% script to retrieve instrument marker by name not last file

%% Defining variables
scriptDir    = fileparts( mfilename('fullpath') );   % e.g. .../projectRoot/scripts

rawDir = fullfile(scriptDir, '..', 'data/raw/Tim_TMS/');  % go up, then into 'functions'
organizedDir = fullfile(scriptDir, '..', 'data/organized/Tim_TMS/');  % go up, then into 'functions'

subjects = 401:441;

%% you may choose what to run on each steps

% step 1: Create new organized Directory for all EntryTarget, InstrumentMarkers and TMSTrigger files
% organizeData(rawDir, organizedDir, subjects);

% step 2: return the filepaths of the organized directory at step 1
size_threshold = 100000000; % TMSTrigger exceeding 100KB are selected

pathTable = extractPaths(organizedDir, subjects, size_threshold);

% ToDo: Instrumentmarker by names
 
% step 3: 
% Use the should-be template (when the session actually take place) to determine 
% the valid relevant file from the list of paths generated in step 2.
% in the following, it's from Google calendar

% finalizing the template here and should-be status here "scheduled",
% "Finished", "Not Scheduled", "Withdrawn".

SessPerSub = 3;
validTable_wide = readtable('valid_session_gcal_wide.csv', Delimiter=","); % u may skip this step if you have long data directly
validTable_long = validSessWide2Long(validTable_wide, SessPerSub); % alternatively: validTable_long = readtable('valid_session_gcal_long.csv', Delimiter=",");

% specifies a two-column table for conditions and their types. This help
% determine what coordinates files would be used for contrast

condMap = table( ...
    ["Vertex"; "iTBS"; "cTBS"], ...                            % Column #1: Condition
    { ["EntryTarget","InstrumentMarkers","TMSTrigger"];        % Column #2: Types
      ["InstrumentMarkers","TMSTrigger"];
      ["InstrumentMarkers","TMSTrigger"] }, ...
    'VariableNames', {'Condition','Types'} ...
);

all_valid_files = mapSessionFiles(validTable_long, pathTable, condMap);

for i=1:height(all_valid_files) % replaces manualReplaceFile function
    if all_valid_files.ID(i) == "CC434" && all_valid_files.Date(i) == 20250506 && all_valid_files.Type(i) == "EntryTarget"
        all_valid_files.FullPath{i} = strcat(organizedDir, "434\Sessions\Manually_created_sessions\EntryTarget\EntryTarget20250909111111111.xml");
        all_valid_files.Status{i} = 'Manual Fill';
    elseif all_valid_files.ID(i) == "CC432" && all_valid_files.Date(i) == 20250527 && all_valid_files.Type(i) == "InstrumentMarkers"
        all_valid_files.FullPath{i} = strcat(organizedDir, "432\Sessions\Manually_created_sessions\InstrumentMarkers\InstrumentMarker20250527145647608.xml");
        all_valid_files.Status{i} = 'Manual Fill';
    elseif all_valid_files.ID(i) == "CC406" && all_valid_files.Date(i) == 20250122 && all_valid_files.Type(i) == "InstrumentMarkers"
        all_valid_files.FullPath{i} = strcat(organizedDir, "406\Sessions\Manually_created_sessions\InstrumentMarkers\InstrumentMarker20250122105915000.xml");
        all_valid_files.Status{i} = 'Manual Fill';
    elseif all_valid_files.ID(i) == "CC408" && all_valid_files.Date(i) == 20250124 && all_valid_files.Type(i) == "EntryTarget"
        all_valid_files.FullPath{i} = strcat(organizedDir, "408\Sessions\Manually_created_sessions\EntryTarget\EntryTarget20250124095420474.xml");
        all_valid_files.Status{i} = 'Manual Fill';
    elseif all_valid_files.ID(i) == "CC415" && all_valid_files.Date(i) == 20250324 && all_valid_files.Type(i) == "InstrumentMarkers"
        all_valid_files.FullPath{i} = strcat(organizedDir, "415\Sessions\Manually_created_sessions\InstrumentMarkers\InstrumentMarker20250324160949000.xml");
        all_valid_files.Status{i} = 'Manual Fill';
    end
end

% step 4: convert from f3 extract coordinate from the valid files

all_files_with_coord = readCoordFromFiles(all_valid_files);

% step 5: 

% since iTBS and cTBS should share the same InstrumentMarker, but sometimes
% only one of them appears, insmarker_equivalence establish such equivalence

all_files_with_coord_eq = insmarker_equivalence(all_files_with_coord, ["iTBS", "cTBS"]); % optional; make a copy of corresponding instrument marker if it doesn't already exist

% every distance condition --------------------------
condSpec = struct( ...
    'CalculatedEntry', struct('RefType','EntryTarget',  'RowCond',"Vertex", 'idx',1 ), ...
    'PresumedTarget',  struct('RefType','EntryTarget',  'RowCond',"Vertex", 'idx',2 ), ... % idx 2 is defaulted for Target locations
    'MarkerEntry',   struct('RefType','InstrumentMarkers','RowCond',"Vertex", 'idx',1 ), ...
    'ImportedMarker_iTBS',    struct('RefType','InstrumentMarkers','RowCond',"iTBS", 'idx',1 ), ...
    'ImportedMarker_cTBS',    struct('RefType','InstrumentMarkers','RowCond',"cTBS", 'idx',1 ));

[TMS_Quality_raw, TMS_Quality_sum] = calculateCoords(all_files_with_coord_eq, condSpec);

fprintf('finished!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TMS_Quality_raw, TMS_Quality_sum] = calculateCoords(all_files_with_coord, condSpec)

    % retrieve the conditions name
    distNames = fieldnames(condSpec);
    
    % split the big table
    
    T = all_files_with_coord;
    isTrigg = T.Type == "TMSTrigger";
    Trigg_T = T(isTrigg,:);
    Ref_T = T(~isTrigg,:);
    
    % Group by subject × date
    
    [G,subj,date] = findgroups( T.ID , T.Date );   % indexing uniques sessions
    nGroups = max(G);
    
    TMS_Quality_sum = table();
    TMS_Quality_raw = table();
    
    % iterate over five conditions
    for g = 1:nGroups % for each unique sessions 
        idxSes = (G==g); % indexes all entries of current unique sessions
        trig = T(idxSes & isTrigg, :); % extract current trigger entry
        ref = T(idxSes & ~isTrigg, :); % extract current non-trigger entry/entries
    
        for h = 1:height(ref) % for vertex, there is two ref. points
            current_ref = ref(h,:);

            for k = 1:numel(distNames) % processing a ref-trg comparison
                nm = distNames{k}; % Ref Type
                spec = condSpec.(nm); % assess Ref Type specs
                
                refRow = table();
                
                % pick the reference row
                refRow = current_ref(current_ref.Type==spec.RefType ... % e.g. EntryTarget
                                   & current_ref.Cond==spec.RowCond, :); % e.g.. Vertex
                if isempty(refRow),  continue, end                 % entry not match to ref → skip
                % pick the trigger row
                trgRow = trig(trig.Cond==spec.RowCond ...
                            & trig.Type=="TMSTrigger", :); % instead of spec.RefType
                if isempty(trgRow),  continue, end                 % entry not match to ref → skip
    
                % -> only the valid rows remain now
    
                % handling missing file for reference row and trigger row
                if trgRow.Date == refRow.Date
                    compDate = trgRow.Date;
                end
                if ~ismember(refRow.Status, ["Finished", "Manual Fill"]) % filter ref coord missing
                    outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
                                         "ref-" + refRow.Status);
                    TMS_Quality_sum = [TMS_Quality_sum;outRow];
                    
                    rawRow = makeRawRow(subj(g), string(nm), NaN, NaN, NaN, NaN, "ref-" + refRow.Status);
                    TMS_Quality_raw = [TMS_Quality_raw; rawRow];
                    continue, 
                end                 % have entry but empty row → skip
                if ~ismember(trgRow.Status, ["Finished", "Manual Fill"]) 
                    outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
                                         "trg-" + trgRow.Status);
                    TMS_Quality_sum = [TMS_Quality_sum;outRow];                
                    
                    rawRow = makeRawRow(subj(g), string(nm), NaN, NaN, NaN, NaN, "ref-" + refRow.Status);
                    TMS_Quality_raw = [TMS_Quality_raw; rawRow];
                    continue                 % have entry but empty row → skip
                end 
                % -> note refRow missing before trgRow
    
                % loading the coordinates 
                refCoord = [refRow.x_axis{1}(spec.idx) ... % inside x_axis content, select the idx's element
                            refRow.y_axis{1}(spec.idx) ...
                            refRow.z_axis{1}(spec.idx)];
                trgCoord = [trgRow.x_axis{:}, trgRow.y_axis{:}, trgRow.z_axis{:}];
    
                % handling RAS coordinate if any
                [tf_refCoord, tf_refCoordSys] = RAS2LPS(refCoord, refRow.coord_system); % convert to same system                            
                [tf_trgCoord, tf_trgCoordSys] = RAS2LPS(trgCoord, trgRow.coord_system); % convert to same system
                
                % skipping MNI files, deciding on output CoordSys name
                CoordCompare = strcat( tf_refCoordSys, "_", tf_trgCoordSys );
    
                % comparing distances, append valid row
                squaredDifferences = (tf_trgCoord - tf_refCoord).^2;

                eucdistances = []; % empty from previous loop

                for h = 1:height(squaredDifferences)
                    eucdistances(h,1) = sqrt(sum(squaredDifferences(h,:))); % row-wise addition
                end
                 
                ed = eucdistances;
    
                outRow = makeQualRow(subj(g), string(nm), compDate, mean(ed), median(ed), max(ed), ...
                                     min(ed), numel(ed), CoordCompare, refRow.Status);   % ref cuz usually ref need manual replacement
                TMS_Quality_sum = [TMS_Quality_sum; outRow];

                rawRow = makeRawRow(subj(g), string(nm), {tf_refCoord}, {tf_trgCoord}, {ed}, numel(ed), refRow.Status); % {} for cell, otherwise can't fit into the table
                TMS_Quality_raw = [TMS_Quality_raw; rawRow];
            end
        end
    end
    disp('finish executing calculateCoords.m\n')
end

% function TMS_Quality_sum = calculateCoords(all_files_with_coord, condSpec)
% 
%     % retrieve the conditions name
%     distNames = fieldnames(condSpec);
% 
%     % split the big table
% 
%     T = all_files_with_coord;
%     isTrigg = T.Type == "TMSTrigger";
%     Trigg_T = T(isTrigg,:);
%     Ref_T = T(~isTrigg,:);
% 
%     % Group by subject × date
% 
%     [G,subj,date] = findgroups( T.ID , T.Date );   % indexing uniques sessions
%     nGroups = max(G);
% 
%     TMS_Quality = table();
% 
%     % iterate over five conditions
%     for g = 1:nGroups % for each unique sessions 
%         idxSes = (G==g); % indexes all entries of current unique sessions
%         trig = T(idxSes & isTrigg, :); % extract current trigger entry
%         ref = T(idxSes & ~isTrigg, :); % extract current non-trigger entry/entries
% 
%         for h = 1:height(ref) % for vertex, there is two ref. points
%             current_ref = ref(h,:);
% 
%             for k = 1:numel(distNames) % processing a ref-trg comparison
%                 nm = distNames{k}; % Ref Type
%                 spec = condSpec.(nm); % assess Ref Type specs
% 
%                 refRow = table();
% 
%                 % pick the reference row
%                 refRow = current_ref(current_ref.Type==spec.RefType ... % e.g. EntryTarget
%                                    & current_ref.Cond==spec.RowCond, :); % e.g.. Vertex
%                 if isempty(refRow),  continue, end                 % entry not match to ref → skip
%                 % pick the trigger row
%                 trgRow = trig(trig.Cond==spec.RowCond ...
%                             & trig.Type=="TMSTrigger", :); % instead of spec.RefType
%                 if isempty(trgRow),  continue, end                 % entry not match to ref → skip
% 
%                 % -> only the valid rows remain now
% 
%                 % handling missing file for reference row and trigger row
%                 if trgRow.Date == refRow.Date
%                     compDate = trgRow.Date;
%                 end
%                 if ~ismember(refRow.Status, ["Finished", "Manual Fill"]) % filter ref coord missing
%                     outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
%                                          "ref-" + refRow.Status);
%                     TMS_Quality = [TMS_Quality;outRow];
%                     continue, 
%                 end                 % have entry but empty row → skip
%                 if ~ismember(trgRow.Status, ["Finished", "Manual Fill"]) 
%                     outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
%                                          "trg-" + trgRow.Status);
%                     TMS_Quality = [TMS_Quality;outRow];                
%                     continue                 % have entry but empty row → skip
%                 end 
%                 % -> note refRow missing before trgRow
% 
%                 % loading the coordinates 
%                 refCoord = [refRow.x_axis{1}(spec.idx) ... % inside x_axis content, select the idx's element
%                             refRow.y_axis{1}(spec.idx) ...
%                             refRow.z_axis{1}(spec.idx)];
%                 trgCoord = [trgRow.x_axis{:}, trgRow.y_axis{:}, trgRow.z_axis{:}];
% 
%                 % handling RAS coordinate if any
%                 [tf_refCoord, tf_refCoordSys] = RAS2LPS(refCoord, refRow.coord_system); % convert to same system                            
%                 [tf_trgCoord, tf_trgCoordSys] = RAS2LPS(trgCoord, trgRow.coord_system); % convert to same system
% 
%                 % skipping MNI files, deciding on output CoordSys name
%                 CoordCompare = strcat( tf_refCoordSys, "_", tf_trgCoordSys );
% 
%                 % comparing distances, append valid row
%                 distances = tf_trgCoord - tf_refCoord; % xyz to xyz
% 
%                 for t=1:length(distances)
%                     euclidean_distances(t,:) = sqrt(sum(distances(t,:).^2));
%                 end
%                 ed = euclidean_distances;
% 
%                 outRow = makeQualRow(subj(g), string(nm), compDate, mean(ed), median(ed), max(ed), ...
%                                      min(ed), numel(ed), CoordCompare, refRow.Status);   % ref cuz usually ref need manual replacement
% 
%                 TMS_Quality = [TMS_Quality;outRow];
%             end
%         end
%     end
%     disp('finish executing calculateCoords.m\n')
% end

function row = makeQualRow(id, cond, compDate, avg, med, maxv, minv, n, coordSysTxt, statusTxt)
    row = table(id, cond, compDate, avg, med, maxv, minv, n, coordSysTxt, statusTxt, ...
        'VariableNames', {'ID', 'distCond', 'compDate', 'avg', 'med', 'max', 'min', 'n', 'CoordSys', 'Status'});
end

function row = makeRawRow(a,b,c,d,e,f,g)
    row = table(a,b,c,d,e,f,g, ...
        'VariableNames', {'ID', 'distCond', 'refCoord', 'trgCoord', 'eucldist', 'n', 'Status'});
end

function [outCoords, xform] = RAS2LPS(coords, coord_system)
    if coord_system == "RAS"
        % flip X and Y (leave Z untouched) convert from RAS to LPS
        outCoords = [-coords(:,1), -coords(:,2),  coords(:,3)];
        xform     = "RAS2LPS";        
    else % for LPS and MNI cases
        xform = coord_system; % initialize, otherwise return nth
        outCoords = coords;% initialize, otherwise return nth
    end
end

function pathTable = extractPaths(organizedDir, subjectIDs, size_threshold)
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
    MIN_FILE_SIZE_THRESHOLD = size_threshold;
    
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
    outDir = fullfile('..','data','processed');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    
    outputFileName = '../data/processed/step_1_all_extracted_path.csv';
    writetable(pathTable, outputFileName);
    disp(['Results saved to ', outputFileName]);
end

function all_files_with_coord_equivalenced = insmarker_equivalence(all_files_with_coord, markers)
    % 1) pull out each sub-table --------------------------------------------
    % markers = string(markers);
    results = struct();
    for mk = markers
        idx = all_files_with_coord.Cond == mk & ...
              all_files_with_coord.Type == "InstrumentMarkers";
        results.(mk) = all_files_with_coord(idx,:);
    end
    
    % 2) compute union -------------------------------------------
    cols = ["x_axis","y_axis","z_axis", "coord_system"];

    nRows = height(results.(markers(1)));
    unionVals = cell(nRows, numel(cols));        % nRows × 3 cell matrix
    
    for j = 1:numel(cols)                        % loop over x/y/z
        c = cols(j);
        colU = cell(nRows,1);                    % starts all empty
    
        for mk = markers                         % scan every marker
            data = results.(mk).(c);
            if ~iscell(data)                    % for coord_system that have string_array
                data = num2cell(data);
            end
                
            needs = cellfun(@isempty, colU);     % rows still unfilled (2nd up)
            colU(needs) = data(needs);           % first non-empty wins
            if ~any(needs), break, end           % done for this column
        end
        unionVals(:,j) = colU;                   % stash
    end
    
    % 3) patch every marker from union ---------------------------------------
    for mk = markers
        T = results.(mk);
        patched = false(nRows,1);
    
        for j = 1:numel(cols)
            c    = cols(j);
            % miss = …                            % you already have this
            uCol = unionVals(:,j);                % column j of the cell union
            if ~iscell(T.(c))                     % this marker holds strings
                uCol = vertcat(uCol{:});          % unwrap back to string array
            end
            
            miss = cellfun(@isempty, T.(c));

            if any(miss)
                T.(c)(miss) = unionVals(miss,j);
                patched     = patched | miss;
            end
        end
    
        % provenance tag (optional but handy)
        if iscell(T.Filename),  T.Filename = string(T.Filename);  end
        T.Filename(patched) = T.Filename(patched) + " FILLED_from_union";
    
        results.(mk) = T;                      % put back
    end
    
    % 4) write everything back to the big table -----------------------------
    out = all_files_with_coord;          
    
    % make the column classes match (string everywhere)
    out.coord_system = string(out.coord_system);
    out.Filename     = string(out.Filename);
    
    for mk = markers
        results.(mk).coord_system = string(results.(mk).coord_system);
        results.(mk).Filename     = string(results.(mk).Filename);
    
        mask = out.Cond == mk & out.Type == "InstrumentMarkers";
        out(mask ,:) = results.(mk);    % now the row-slice assignment works
    end
    
    all_files_with_coord_equivalenced = out;
    disp('finish executing insmarker_equivalence.m\n')
end

function all_valid_files = mapSessionFiles(validSessionLong, pathTable, condMap)
    % to make things simple, we could turn the subj id and date into long
    % format, then use the join function to map it with the relevant file path,
    % extracting the path info only
    
    % conditions csv + path csv -> all valid files mat
    % important % 
    % ########to be added , initialized all possible combinations of files, 
    % before commiting to join, thereby can list out reason more easily
    % --> adding status column
    % #######################
    
        % condMap is a table with columns:
        %   condMap.Condition  – e.g. "Vertex","iTBS","cTBS"
        %   condMap.Types      – a cell array; each cell is a string array of Types
    
    % idx array for duplication
    validSessionLong.Cond = string(validSessionLong.Cond); % turn into string or categorical
    
    type_row = [];
    % initialize FileTypes per session Type
    for r=1:height(condMap)
        currCond = condMap.Condition(r); % e.g. vertex
        currTypeList = condMap.Types{r}; % e.g. instrumentmarker, entrytarget, TMStrigger

        current_idx_array = currCond == validSessionLong.Cond;
    
        base = validSessionLong(current_idx_array, :); % condition specific duplicates

        count = height(base); % count of session per condition

        dupIdx = repelem((1:count), length(currTypeList)); % design the numbering for the new rows location
        type_rows = base(dupIdx,:); % turn number array to rows
    
        type_rows.Type = repmat(currTypeList, 1, count)';
        type_row = [type_row; type_rows];
    end
        
    % to do: join rows with links, and then if no files, return missing in status.
    % to do
    
    file_info = pathTable;
    
    % prepare for join - make value comparable
    keyVars = {'ID' 'Date', 'Type'};
    file_info.ID = file_info.Subject;
    file_info.Type = string(file_info.Type);
    type_row.Date = type_row.TMS_Date;
    file_info.Date = str2double(file_info.Date);
    file_info.Subject = [];
    
    % S_Long.Date = cell2mat(S_Long.Date);
    % S_Long.Date = cellfun(cell2mat, S_Long.Date); converting double to cell
    % array very hard that it's stupid
    
    % join for filename
    S_Long = outerjoin(type_row, file_info, ...
                        'Type', 'left', ...
                        'key', keyVars, ...
                        'MergeKeys', true);
    
    % 'MergeKeys', true) file_info(:, [keyVars 'Filename']), ...
    
    mask_missing = S_Long.Filename == "" & S_Long.Status == "Finished";
    
    S_Long.Status(mask_missing) = "Missing File";
    
    all_valid_files_unfiltered = S_Long;
    
    % remove last redundent file, 
    mask = contains(S_Long.Type,"EntryTarget") & S_Long.Cond~="Vertex";
    
    % drop those rows:
    S_Long(mask, :) = [];
    
    all_valid_files = S_Long;
    
    % initialize Cond_Type
    all_valid_files.Cond_Type = strcat(all_valid_files.Cond, '_', all_valid_files.Type);
    
    
    % (2a) Using ismember → rows in A whose Filename is NOT in B:
    mask = ~ismember(all_valid_files_unfiltered.Filename, all_valid_files.Filename);
    onlyInA = all_valid_files_unfiltered(mask,:);
    
    % remove redudent EntryTarge
    save("all_valid_files.mat", 'all_valid_files')
    disp('finish executing mapSessionFiles.m\n')

end

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

function all_files_with_coord = readCoordFromFiles(all_valid_files)

% reading x y z for all files, and appending it to the table would have
% been more straightforward and organized.

% since we are not always using the default 1x1 array (but n x 1 array) we need to first
% declare it 
    all_valid_files.x_axis = cell(height(all_valid_files),1);
    all_valid_files.y_axis = cell(height(all_valid_files),1);
    all_valid_files.z_axis = cell(height(all_valid_files),1);
    
    
    for f = 1:height(all_valid_files)
        current_sub = all_valid_files.ID(f);
        
        current_path = string(all_valid_files.FullPath(f)); %path is cell originally, need turn to string for readstruct
        
        if current_path == ""
            continue
        end
    
        info = readstruct(current_path); 
       
        current_Type = all_valid_files.Type(f);
    
        if strcmp(current_Type, 'InstrumentMarkers')
            all_valid_files.coord_system(f) = info.coordinateSpaceAttribute;
            current_Matrix4D = info.InstrumentMarker.Marker.Matrix4D; % helper
            all_valid_files.x_axis{f} = current_Matrix4D.data03Attribute;
            all_valid_files.y_axis{f} = current_Matrix4D.data13Attribute;
            all_valid_files.z_axis{f} = current_Matrix4D.data23Attribute;
        
        elseif strcmp(current_Type, 'EntryTarget')
            % for i = 1:length(selectedAttribute)
            %     trueFalse = selectedAttribute(i);
            %     if strcmp(trueFalse, "true")
            %         continue % break the loop and keep the current i
            %     end
            % end
            all_sites=extractfield(info.Entry, 'selectedAttribute');
        
            % extract stim entry
            stim_entry=find(ismember([all_sites{:}],'true')==1); % info.entry click --> marker--> have name for either     
            if length(stim_entry)>1
                error('More than one stimulation sites')
            end        
            % extract stim target
            stim_target=find(ismember([all_sites{:}],'true')==1); % info.entry click --> marker--> have name for either     
            if length(stim_target)>1
                error('More than one stimulation sites')
            end
    
            % store Entry, Target  in cell array
            all_valid_files.coord_system(f) = info.coordinateSpaceAttribute;
            %     coord_system = extractfield(info, 'coordinateSpaceAttribute'); % E.g., 'RAS', 'LAS'
    
    
            all_valid_files.x_axis{f} = [info.Entry(stim_entry).Marker.ColVec3D.data0Attribute; ...
                                            info.Target(stim_target).Marker.ColVec3D.data0Attribute]; 
            all_valid_files.y_axis{f} = [info.Entry(stim_entry).Marker.ColVec3D.data1Attribute; ...
                                            info.Target(stim_target).Marker.ColVec3D.data1Attribute];  
            all_valid_files.z_axis{f} = [info.Entry(stim_entry).Marker.ColVec3D.data2Attribute; ...
                                            info.Target(stim_target).Marker.ColVec3D.data2Attribute];  
        elseif strcmp(current_Type, 'TMSTrigger')
            % store 200 points in cell array
            
            all_valid_files.coord_system(f) = info.coordinateSpaceAttribute;
    
            if ~isfield(info, 'TriggerMarker')
                all_valid_files.coord_system = [];
            else
                all_sites=extractfield(info, 'TriggerMarker');
                numRows = size(all_sites{1}, 2);
    
                if numRows < 200
                    disp(['no enough markers for ', current_ID]);
                end
                
                for k=1:numRows
                    % % Extract the     amplitude
                    % fixed to correctly identify the vertex location [with vertex name]
            
                    % % Extract the coordinate
                    % all_valid_files(k).valueA = all_sites{1, 1}(k).ResponseValues.Value(1).responseAttribute;
                    % all_valid_files(k).amplitudeA = all_sites{1, 1}(k).ResponseValues.Value(6).responseAttribute;
                    x_trig = all_sites{1, 1}(k).Matrix4D.data03Attribute;
                    y_trig = all_sites{1, 1}(k).Matrix4D.data13Attribute;
                    z_trig = all_sites{1, 1}(k).Matrix4D.data23Attribute;
                    
                    % filter out empty coordinates
                    keep_row = ~(x_trig == 0 & y_trig == 0 & z_trig == 0);
    
                    if keep_row 
                        all_valid_files.x_axis{f} = [all_valid_files.x_axis{f}; x_trig]; % append x
                        all_valid_files.y_axis{f} = [all_valid_files.y_axis{f}; y_trig];
                        all_valid_files.z_axis{f} = [all_valid_files.z_axis{f}; z_trig];
                    end
                end
            end
        else
            sprintf('Do not have relevant Type for %s', current_path);
        end
    end
    
    % drop excess column
    all_valid_files(:, 'TMS_Date') = [];
        
    all_files_with_coord = all_valid_files;
    
    save("all_files_with_coord.mat", "all_files_with_coord")
    disp('finish executing readCoordFromFiles.m\n')

end

function valid_session_long = validSessWide2Long(validTable_wide, SessPerSub)
    
    dateVars = arrayfun(@(i) sprintf('TMS%d_Date',  i), 1:SessPerSub, 'UniformOutput', false);
    condVars = arrayfun(@(i) sprintf('TMS%d_Cond',  i), 1:SessPerSub, 'UniformOutput', false);
    
    subj_cond_long = stack( ...
        validTable_wide, dateVars, ...
        'NewDataVariableName',    'TMS_Date',      ...
        'IndexVariableName',      'TMS_Date_Indicator' ...
    );
    
    % created new variable TMS_Date_Indicator.
    subj_cond_long.TMS_Date_Indicator = string(subj_cond_long.TMS_Date_Indicator);
    subj_cond_long.Date_order = cellfun(@(x)replace(x, 'Date', 'Cond'), ...
                                            subj_cond_long.TMS_Date_Indicator, ...
                                            'UniformOutput', false);
    
    idx = (1:height(subj_cond_long))';
    subj_cond_long.Cond = arrayfun( ...
        @(i) subj_cond_long{i, subj_cond_long.Date_order{i}}{1}, ... % T(rowidx, colidx)
        idx, ...
        'UniformOutput', false);
    
    % tag missing date, remove redundent columns
    subj_cond_long.Status = repmat("Finished", height(subj_cond_long), 1);
    
    mask1 = ismissing(subj_cond_long.Cond);  
    subj_cond_long.Status(mask1) = "Withdrawn";
    
    mask2 = isnan(subj_cond_long.TMS_Date) & ~ismissing(subj_cond_long.Cond);  
    subj_cond_long.Status(mask2) = "Not scheduled";
    
    today_dn = datetime('now');
    str_dates = string(subj_cond_long.TMS_Date); % double 2 str
    mask3 = datetime(str_dates, 'InputFormat', 'yyyyMMdd') >= today_dn; % str 2 datetime
    subj_cond_long.Status(mask3) = "Scheduled";
    
    toRemove = [condVars, {'TMS_Date_Indicator'}];
    valid_session_long = removevars(subj_cond_long, toRemove);

    %% CREATE OUTPUT TABLE
    outputFileName = 'valid_session_gcal_long.csv';
    writetable(valid_session_long, outputFileName);
    disp(['Results saved to ', outputFileName]);
    disp('finish executing wideSessWide2Long.m\n')
end
