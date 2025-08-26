function [CsvCompareCoords, CsvSummaryTable] = CsvCompareCoords(all_files_with_coord_eq, csvPath)
%COMPARECOORDS Compare native coords from CSV to coords in all_files_with_coord_eq.
%   CsvCompareCoords = CsvCompareCoords(all_files_with_coord_eq)
%   CsvCompareCoords = CsvCompareCoords(all_files_with_coord_eq, 'native_coords_summary.csv')
%
% Output columns:
%   ID (double)
%   Type (string)                % 'InstrumentMarkers' or 'TMSTrigger'
%   Cond (string)                % copied from all_files_with_coord_eq
%   Date (datetime or string)    % copied from all_files_with_coord_eq
%   Coords_CSV (cell -> 1x3 double)
%   Coords_Table (cell -> Nx3 double)
%   Mean_Distance (double)            % mean Euclidean distance to CSV coord
%   Distances_All (cell -> Nx1 double)

    if nargin < 2 || isempty(csvPath)
        csvPath = 'native_coords_summary.csv';
    end

    % ---------- 1) Read CSV and grab ID + Native coords ----------
    csvT = readtable(csvPath);
    idCol = 'Subject_ID';
    xCol  = 'LPS_X';
    yCol  = 'LPS_Y';
    zCol  = 'LPS_Z';

   % Build a map: ID -> [X Y Z]
    id2coord = containers.Map('KeyType','double','ValueType','any');
    for i = 1:height(csvT)
        thisID = str2double(regexprep(string(csvT.(idCol)(i)), '\D', ''));
        if ~isnan(thisID)
            coords = [csvT.(xCol)(i), csvT.(yCol)(i), csvT.(zCol)(i)];
            id2coord(thisID) = coords;
        end
    end

    %% 2) Prepare collectors
    outID = {}; outType = {}; outCond = {}; outDate = {};
    outCSV = {}; outTable = {}; outDist = []; outDistsAll = {};

       %% 3) Loop through all_files_with_coord_eq
    for r = 1:height(all_files_with_coord_eq)
        % Clean ID (strip CC prefix etc.)
        rawID = string(all_files_with_coord_eq.ID(r));
        numID = str2double(regexprep(rawID, '\D', ''));
        if isnan(numID) || ~isKey(id2coord, numID)
            continue;
        end

        % Only keep relevant types
        thisType = string(all_files_with_coord_eq.Type(r));
        if ~(thisType == "InstrumentMarkers" || thisType == "TMSTrigger")
            continue;
        end

        % Grab CSV coords
        csvCoord = id2coord(numID);

        % Extract coords from this row
        coordsTbl = stackXYZ(all_files_with_coord_eq.x_axis(r), ...
                             all_files_with_coord_eq.y_axis(r), ...
                             all_files_with_coord_eq.z_axis(r));

        % Distances
        if isempty(coordsTbl)
            dists = NaN; dMean = NaN;
        else
            diffs = coordsTbl - csvCoord;
            dists = sqrt(sum(diffs.^2, 2));
            dMean = mean(dists, 'omitnan');
        end

        % Cond / Date
        thisCond = string(all_files_with_coord_eq.Cond(r));
        thisDate = string(all_files_with_coord_eq.Date(r));

        % Append
        outID{end+1,1} = numID;
        outType{end+1,1} = thisType;
        outCond{end+1,1} = thisCond;
        outDate{end+1,1} = thisDate;
        outCSV{end+1,1} = csvCoord;
        outTable{end+1,1} = coordsTbl;
        outDist(end+1,1) = dMean;
        outDistsAll{end+1,1} = dists;
    end

    %% 4) Final table
    CsvCompareCoords = table( ...
        cell2mat(outID), string(outType), string(outCond), string(outDate), ...
        outCSV, outTable, outDist, outDistsAll, ...
        'VariableNames', ...
        {'ID','Type','Cond','Date','Coords_CSV','Coords_Table','Mean_Distance','Distances_All'} ...
    );

    % ---------- 5) Build summary table ----------
    % Remove 'Vertex' rows
    filtered = CsvCompareCoords(CsvCompareCoords.Cond ~= "Vertex", :);

    uniqueIDs = unique(filtered.ID);
    subID = [];
    trigAvg = [];
    markerAvg = [];

    for i = 1:numel(uniqueIDs)
        uid = uniqueIDs(i);
        rows = filtered(filtered.ID == uid, :);

        % TMSTrigger mean of 2 rows
        trigRows = rows(rows.Type == "TMSTrigger", :);
        if height(trigRows) ~= 2
            warning('Subject %d: Expected 2 TMSTrigger rows, found %d', uid, height(trigRows));
        end
        trigMean = mean(trigRows.Mean_Distance, 'omitnan');

        % InstrumentMarkers mean of 2 rows
        markerRows = rows(rows.Type == "InstrumentMarkers", :);
        if height(markerRows) ~= 2
            warning('Subject %d: Expected 2 InstrumentMarkers rows, found %d', uid, height(markerRows));
        end
        markerVals = markerRows.Mean_Distance;
        if height(markerRows) == 2 && abs(markerVals(1) - markerVals(2)) > 1e-6
            warning('Subject %d: InstrumentMarkers values do not match: %f vs %f', uid, markerVals(1), markerVals(2));
        end
        markerMean = mean(markerVals, 'omitnan');

        subID(end+1,1) = uid;
        trigAvg(end+1,1) = trigMean;
        markerAvg(end+1,1) = markerMean;
    end

    CsvSummaryTable = table(subID, trigAvg, markerAvg, ...
        'VariableNames', {'ID','TMSTrigger_Avg','InstrumentMarkers_Avg'});


    fprintf('finished executing CsvCompareCoords.m\n');
end

%% ========== Helper ==========
function XYZ = stackXYZ(xv, yv, zv)
% Normalize x/y/z values from one table row into an Nx3 numeric matrix.
% Handles scalars, vectors, and cell-wrapped numeric arrays.

    x = unwrapNumeric(xv); y = unwrapNumeric(yv); z = unwrapNumeric(zv);

    % Make column vectors
    x = x(:); y = y(:); z = z(:);

    % Handle length mismatches defensively (use the common min length)
    n = min([numel(x), numel(y), numel(z)]);
    if n == 0
        XYZ = [];
        return;
    end

    XYZ = [x(1:n), y(1:n), z(1:n)];
end

function a = unwrapNumeric(v)
% If v is a cell, take its contents; otherwise return numeric as-is.
% If empty or non-numeric, return [].
    if iscell(v)
        if isempty(v)
            a = [];
            return;
        end
        v = v{1};
    end
    if isnumeric(v)
        a = v;
    else
        a = [];
    end
end
