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
    % save("all_valid_files.mat", 'all_valid_files')
    disp('finish executing mapSessionFiles.m\n')

end
