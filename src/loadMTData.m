function neurons = loadMTData(N)
%LOADMTDATA Loads N random .mat files in the MT_Data folder. If N is not
% specified, loads up all of them. If N is a vector of indices, loads up
% those indices to allow for repeatable runs (if any index is too big,
% returns an error).
%   Returns:
%     neurons: array of structures
    if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
        MT_DATA_LOC = 'MT_data';
    else
        MT_DATA_LOC = strcat('..', filesep, 'MT_data');
    end

    % Get all files in MT_data folder
    dircon = dir(MT_DATA_LOC);    % get all files in folder
    files = cell([length(dircon), 1]);
    for i = 1:length(dircon)
        files{i} = dircon(i).name;
    end
    
    % Get the indices of files in the folder that are .mat data files for
    % neuron activity. Files are of the format: "cell_xxxxx.mat"
    inds = contains(files, '.mat') & contains(files, 'cell');
    
    % If N > number of relevant files, produce error
    n_rel = sum(inds);
    if (nargin == 0)
        N = n_rel;
    elseif (any(N > n_rel))
        error('Only %d MT neuron data files are available.', n_rel);
    end
    
    % Take only the relevant files
    files = files(inds);
    
    if (length(N) == 1)
        % Generate N random indices
        fnum = randperm(n_rel, N);
    else
        fnum = N;
        N = length(fnum);
    end
    
    % Load chosen files into neurons, an array of structs
    loads = files(fnum);
    neurons = [];  % Can you pre-initialize a struct array??
    for f = 1:N
        n = load(strcat(MT_DATA_LOC, filesep, cell2mat(loads(f,:))));
        n.data(n.data(:) ~= 0) = 1;  % Make data more binary
        n.name = loads(f,:);
        neurons = [neurons, n]; %#ok<AGROW>
    end
end

