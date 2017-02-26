function neuron = genNeuron(ID, ndirs, stim, reps)
%GENNEURON Generates a single neuron, represented by a struct with fields:
%   experiment: Neuron ID, as a passed-in string
%   dirs: Vector of stimulus directions, spread evenly over ndirs steps
%     from -90 to 90.
%   maxrep: ???
%   data: Using the stim vector (representing the stimulus), generates
%     reps repetitions of firing data, represented by a matrix of
%     dimensions: length(stim) x dirs+1 x reps
%   tdata: matrix of spike times
    T = length(stim);

    % Get stimulus directions vector. Use ndirs-1 since the direction
    % angles themselves describe the "edges" of the range of directions
    dirs = -90 : 180 / (ndirs - 1) : 90;
    
    % Maxrep???
    mr = reps;
    
    % Generate neuron firing data
    data = zeros(T, dirs+1, reps);
    
    % Spike times
    tdata = zeros(size(data));
    
    % Build neuron struct for output
    neuron = struct('experiment', ID, 'dirs', dirs, 'maxrep', mr, ...
                    'data', data, 'tdata', tdata);
end

