function all_valid_files = manualReplaceFile(all_valid_files, overwrite_csv)
%MANUALREPLACEFILE Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%   add back the feature for specifying removal file TYPE, to account for indexing in the empty cases
    if height(overwrite_csv) < 1
        return
    end
    config = pipeline_config();

    for r =1:height(overwrite_csv)

        cur_subj = overwrite_csv.Subject{r};
        cur_date = overwrite_csv.Date(r);    

        if ~isnan(overwrite_csv.FullPath(r)) % if not empty/NaN, do search
            cur_fpath = [config.dirs.organized, '\', overwrite_csv.FullPath{r}];
    
            [~, baseName, ext] = fileparts(cur_fpath);
            cur_fname = strcat(baseName, ext);
        
            % check date matches
            dateStrs = regexp(cur_fname,'\d{8}','match','once');
            fileDate = str2double(dateStrs);
            if cur_date ~=  fileDate
                fprintf('warning: mismatch input date %d and file date %d]\n', cur_date, fileDate)
            end
    
            type = regexp(cur_fname, '^[^\d]*', 'match', 'once');
    
            if strcmp(type, 'InstrumentMarker')
                type = 'InstrumentMarkers'; % Localite internal naming excentricities (add s)
            end
            targetRow = strcmp(all_valid_files.ID, cur_subj) & all_valid_files.Date == cur_date & strcmp(all_valid_files.Type, type);
            if ~any(targetRow)
                fprintf('critical warning: target row not found\n')
            end
            all_valid_files.FullPath(targetRow) = {cur_fpath};
            all_valid_files.Status(targetRow) = "Manual Fill";
    
            all_valid_files.Filename(targetRow) = {cur_fname}; % concated chars and put into a cell  

        else % cur_path empty, remove the existing file reference to nothing

            fprintf("removing file reference for %d on %d\n", cur_subj, cur_date)


            targetRow = strcmp(all_valid_files.ID, cur_subj) & all_valid_files.Date == cur_date & strcmp(all_valid_files.Type, type);
            if ~any(targetRow)
                fprintf('critical warning: target row not found\n')
            end
            all_valid_files.FullPath(targetRow) = '';
            all_valid_files.Status(targetRow) = "Manual Removed";
    
            all_valid_files.Filename(targetRow) = ''; % concated chars and put into a cell  
        end
    end
    return
end

