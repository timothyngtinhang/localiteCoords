function checkDurationMatch(all_files, config)
%CHECKDURATIONMATCH Checks if durations match expected values per condition.
%   This function compares the duration of each TMS session against the
%   expected duration and tolerance values defined in the config file.
%
% Inputs:
%   all_files - Table containing at least these columns:
%               'Cond'         : condition name (string or categorical)
%               'duration_sec' : duration in seconds (numeric)
%
% Behaviour:
%   Prints warnings in the command window for:
%     - Unknown conditions
%     - Durations outside expected ranges (with ± tolerance)
%
% Notes:
%   - Skips rows with missing (NaN) or zero duration.
%   - Expected durations:
%       cTBS and Vertex: ~40s ± 5s
%       iTBS: ~188s ± 5s

    % Get the expectedDurations struct from config
    expectedDurations = config.durationCheck.expectedDurations;
    hasSuspicious = false; %NEW: Track if any suspicious duration was found

    for idx = 1:height(all_files)
        cond = string(all_files.Cond(idx));
        duration = all_files.duration_sec(idx);

        % Skip if duration is missing or zero
        if isnan(duration) || duration == 0
            continue
        end

%         % Determine expected duration and tolerance
%         switch cond
%             case {"cTBS", "Vertex"}
%                 expected = 40;
%                 tolerance = 1;
%             case "iTBS"
%                 expected = 188;
%                 tolerance = 1;
%             otherwise
%                 fprintf('Warning: Unknown condition "%s" in row %d\n', cond, idx);
%                 continue
%         end

        % Check if the condition exists in config
        if isfield(expectedDurations, cond)
            expected = expectedDurations.(cond).expected;
            tolerance = expectedDurations.(cond).tolerance;
        else
            fprintf('Warning: Unknown condition "%s" in row %d\n', cond, idx);
            continue
        end

        % Check if duration is within expected range
        if duration < (expected - tolerance) || duration > (expected + tolerance)
            fprintf('[Duration Check] Condition "%s" in row %d has suspicious duration: %.1f sec (expected around %d sec)\n', ...
                cond, idx, duration, expected);
            hasSuspicious = true; % Flag it
        end
    end

    % NEW: Report if all durations are fine
    if ~hasSuspicious
        fprintf('[Duration Check] All durations are within the expected range.\n');
    end
end