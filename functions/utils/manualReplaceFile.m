function all_valid_files = manualReplaceFile(all_valid_files, overwrite_csv)
% MANUALREPLACEFILE Replaces or clears file paths manually based on user-provided table.
%
% This function processes the `overwrite_csv` table which may contain:
% - Valid new file paths to manually overwrite existing ones
% - Empty paths that signal removal of invalid entries
%
% It updates the corresponding entries in `all_valid_files` based on Subject, Date, and Type.
    
    if height(overwrite_csv) < 1
        return
    end

     % --- Normalize FullPath column to ensure safe indexing ---
    if iscell(overwrite_csv.FullPath)
        % Already in correct format: cell array of char vectors
    elseif isstring(overwrite_csv.FullPath)
        % Convert string array to cell array of char vectors
        overwrite_csv.FullPath = cellstr(overwrite_csv.FullPath);
    elseif isnumeric(overwrite_csv.FullPath)
        % If it's numeric (e.g. full of NaNs), replace with empty strings
        overwrite_csv.FullPath = repmat({''}, height(overwrite_csv), 1);
    else
        % Unexpected format → throw error
        error('Unsupported format for FullPath column');
    end

    config = pipeline_config();

    for r = 1:height(overwrite_csv)

        cur_subj = overwrite_csv.Subject{r};
        cur_date = overwrite_csv.Date(r);
        cur_path_raw = overwrite_csv.FullPath{r}; % Raw path provided (may be empty)
        
        %NEW:Default to empty filename and path
        cur_fname = "";
        cur_fpath = "";


      %If path is provided, build full path and extract filenametype
        if ~isempty(cur_path_raw)
            cur_fpath = fullfile(config.dirs.organized, cur_path_raw); % Combine base directory with relative path to get full path
            [~, baseName, ext] = fileparts(cur_fpath); % Extract base filename and extension
            cur_fname = strcat(baseName, ext); % Combine to get full filename

            % Extract date from filename and compare to provided date
            dateStrs = regexp(cur_fname,'\d{8}','match','once');
            fileDate = str2double(dateStrs);
            if cur_date ~=  fileDate
            fprintf('warning: mismatch input date %d and file date %d]\n', cur_date, fileDate)
            end
        
            % Extract type (e.g., "EntryTarget") from filename prefix
            type = regexp(cur_fname, '^[^\d]*', 'match', 'once');
        else
            %If path is empty, still try to use type from csv directly
            type = overwrite_csv.Type{r};
        end

        if strcmp(type, 'InstrumentMarker')
            type = 'InstrumentMarkers'; % Localite internal naming excentricities (add s)
        end
        
       % Identify the target row(s) in `all_valid_files` that match subject, date, and type
       targetRow = strcmp(all_valid_files.ID, cur_subj) & ...
                    all_valid_files.Date == cur_date & ...
                    strcmp(all_valid_files.Type, type);

        if ~any(targetRow)
            fprintf('critical warning: target row not found for %s %d %s \n', cur_subj, cur_date, type);
            continue
        end

        % === Case 1: Path was blank → remove file reference ===
        if isempty(cur_path_raw)
            all_valid_files.FullPath(targetRow) = {''};
            all_valid_files.Filename(targetRow) = {''};
            all_valid_files.Status(targetRow) = "Manual Fill (missing)";
        
        % === Case 2: Path was provided → manually set file ===
        else
            all_valid_files.FullPath(targetRow) = {cur_fpath};
            all_valid_files.Filename(targetRow) = {cur_fname};
            all_valid_files.Status(targetRow) = "Manual Fill"; % concated chars and put into a cell  
        end
    end
end
