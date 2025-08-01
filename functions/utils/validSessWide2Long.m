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
