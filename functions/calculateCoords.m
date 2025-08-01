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
                    
                    rawRow = makeRawRow(subj(g), string(nm), {NaN}, {NaN}, {NaN}, NaN, "ref-" + refRow.Status);
                    TMS_Quality_raw = [TMS_Quality_raw; rawRow];
                    continue, 
                end                 % have entry but empty row → skip
                if ~ismember(trgRow.Status, ["Finished", "Manual Fill"]) 
                    outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
                                         "trg-" + trgRow.Status);
                    TMS_Quality_sum = [TMS_Quality_sum;outRow];                
                    
                    rawRow = makeRawRow(subj(g), string(nm), {NaN}, {NaN}, {NaN}, NaN, "ref-" + refRow.Status);
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
                
                if norm(tf_refCoord) < 0.1 % basically cases when refCoord = [0, 0, 0] or very closed to it / L2norm
                    outRow = makeQualRow(subj(g), string(nm), compDate, NaN, NaN, NaN, NaN, NaN, NaN, ...
                                         "refCoord ~ [0,0,0]");
                    TMS_Quality_sum = [TMS_Quality_sum;outRow];                

                    rawRow = makeRawRow(subj(g), string(nm), {NaN}, {NaN}, {NaN}, NaN, "refCoord = [0,0,0]");
                    TMS_Quality_raw = [TMS_Quality_raw; rawRow];
                    continue
                end

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
