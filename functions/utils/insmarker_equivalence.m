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

