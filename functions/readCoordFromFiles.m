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
            % sometimes instrumentmarker contains two coordinate, choose
            % the one with 'selected' set to 'true'
            % but sometimes instrumentmarker file have only 1 coordinate,
            % with 'selected' set to 'false', so selectedMask is used only
            % when there is more than 1 coordinate
            if length(info.InstrumentMarker) > 1
                selectedMask = strcmp([info.InstrumentMarker.selectedAttribute], "true");
                current_Matrix4D = info.InstrumentMarker(selectedMask).Marker.Matrix4D;
            else
                current_Matrix4D = info.InstrumentMarker.Marker.Matrix4D;
            end

            % current_Matrix4D = info.InstrumentMarker.Marker.Matrix4D; % helper

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
                    disp(['no enough markers for ', current_sub]);
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
    
    % save("all_files_with_coord.mat", "all_files_with_coord")
    disp('finish executing readCoordFromFiles.m\n')

end

% conditions = {'Vertex', 'iTBS', 'cTBS'};
% types = {'EntryTarget', 'InstrumentMarkers', 'TMSTrigger'};
% 
% dataCoordSystem = struct('ID', {}, 'coord_system', {}, 'x_axis', {}, 'y_axis', {}, 'z_axis', {});
% 
% for c = 1:length(conditions)
%     current_cond = conditions(c);
%     for t = 1:length(types)
%         current_type = types(t);
% 
%         if strcmp(current_type, 'InstrumentMarkers')
%             for s = 1:height(all_valid_files.ID)
%                 current_ID = all_valid_files.ID(s);
%             end
%         end
%     end
%     % no need remove missing files since source table already used inner
%     % join to remove missing
% end







