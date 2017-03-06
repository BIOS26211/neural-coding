function code = getCoding(neurons)
%GETCODING For the given input neurons, returns the neural coding for that
% population in the structure code. Code contains the following fields:
%   tbins: number of time bins (for all neurons)
%   reps: number of reps used. Equal to the number of reps used for the
%     individual neuron with the fewest reps.
%   rnum: repetition number
%   code: a 4D matrix of neural codings. The first two dimensions represent
%     individual neurons (rows) and time bins (columns), respectively.
%     Dimension 3 is for the 14 different stimulus directions, while the
%     last is for the different trials. Draws a number of trials equal to
%     the least number of trials over all the neurons presented.
%   length: number of neurons in each codeword
%   size: number of unique codewords in the code
%   words: all unique words
%   weights: Hamming weights of each word (number of 1s in the word)
%   sparsity: average weight of all codewords in the code
%   redundancy: quantifies the idea that typically more neurons are used
%     than would be necessary to encode a given set of stimuli. This can be
%     interpreted as redundancy*length = number of redundant neurons in the
%     code.

    N = length(neurons);  % number of neurons
    if (N == 0)
        error('Neural coding requires at least one neuron.')
    end
    
    % Get parameters
    tbins = length(neurons(1).data);  % Number of time bins
    reps = min([neurons(:).maxrep]);    % Number of reps
    
    % Set to structure
    code.tbins = tbins;
    code.reps = reps;
    code.rnum = zeros(N, 14, reps);
    code.code = zeros(N, tbins, 14, reps);
    code.length = N;
    
    % Hard code directions to make sure we select these
    dirs = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90];
    
    % Build neural code
    % Iterate over each neuron
    for n = 1:N
        % Iterate over stimulus direction
        [dInd, cd] = ismember(neurons(n).dirs, dirs);
        for d = 1:length(dInd)
            if (dInd(d) == 0)
                continue;
            end
            
            rnums = randperm(neurons(n).maxrep, reps);
            code.rnum(n, cd(d), :) = rnums;
            % Iterate over repetitions
            for r = 1:reps
                c = neurons(n).data(:, d, rnums(r));
                c(c(:) ~= 0) = 1;  % Make sure non-zero values = 1
                code.code(n, :, cd(d), r) = c;
            end
        end
    end
    
    % Unique words
    unis = zeros(tbins * 14 * reps, N);
    u = 1;
    
    % Get unique words
    
    for t = 1:tbins
        for d = 1:length(dirs)
            for r = 1:reps
                unis(u,:) = code.code(:, t, d, r);
                u = u + 1;
            end
        end
    end
   
    code.words = unique(unis, 'rows')';
    code.size = length(code.words);
    code.weights = sum(code.words);
    code.sparsity = sum(code.weights) / (code.size * N);
    code.redundancy = 1 - log2(code.size) / N;
end

