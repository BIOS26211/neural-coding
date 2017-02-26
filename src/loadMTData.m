function neurons = loadMTData()
%LOADMTDATA Loads each .mat file in the MT_Data folder.
%   Returns:
%     neurons: array of structures
    if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
        MT_DATA_LOC = 'MT_data';
    else
        MT_DATA_LOC = strcat('..', filesep, 'MT_data');
    end

    % Get all files in MT_data folder
    files = ls(MT_DATA_LOC);
    
    neurons = [];
    for f = 1:size(files,1)
        if (~isempty(strfind(files(f,:), '.mat')) && ~isempty(strfind(files(f,:), 'cell')))
            % Load and add
            n = load(strcat(MT_DATA_LOC, filesep, files(f,:)));
            neurons = [neurons, n]; %#ok<AGROW>
        end
    end
end

