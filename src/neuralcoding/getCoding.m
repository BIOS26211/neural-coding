function code = getCoding(neurons)
%GETCODING For the given input neurons, returns the neural coding for that
% population in the structure code. Input parameter neurons should be an
% array of structs, with each struct containing information about a single
% neuron.
% Code contains the following fields:
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
%   wordprobs: probability of each word occuring over all time steps and
%     directions.
%   weights: Hamming weights of each word (number of 1s in the word)
%   sparsity: average weight of all codewords in the code
%   redundancy: quantifies the idea that typically more neurons are used
%     than would be necessary to encode a given set of stimuli. This can be
%     interpreted as redundancy*length = number of redundant neurons in the
%     code.

    N = length(neurons);  % number of neurons
    if (N == 0)
        error('Neural coding requires at least one neuron.')
    elseif (~any(strcmp('data', fieldnames(neurons(1)))))
        error('Couldnt get neuron data');
    end
    
    % Get parameters
    tbins = length(neurons(1).data);  % Number of time bins
    reps = min([neurons(:).maxrep]);    % Number of reps
    
    % Set to structure
    code.tbins = tbins;
    code.reps = reps;
    code.length = N;
    
    % Hard code directions to make sure we select these
    dirs = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90];
    code.dirs = dirs;
    nDirs = length(dirs) + 1;
    
    code.rnum = zeros(N, nDirs, reps);
    code.code = zeros(N, tbins, nDirs, reps);
    
    % Build neural code
    % Iterate over each neuron
    for n = 1:N
        % Iterate over stimulus direction
        % +1 to account for null dir
        dInd = ones(1, length(neurons(n).dirs) + 1);
        cd = ones(1, length(neurons(n).dirs) + 1);
        [dInd(2:end), cd(2:end)] = ismember(neurons(n).dirs, dirs);
        % Increments elements in cd except the first and those == 0 by 1
        cd([0, 2:end] & cd > 0) = cd([0, 2:end] & cd > 0) + 1;
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
    unis = reshape(code.code, [N, tbins * nDirs * reps]);
    code.words = unique(unis', 'rows')';
    
    % Get word probabilities
    wordcount = numel(code.code(1, :, 1, :));
    tWC = numel(code.code(1, 1, 1, :));
    code.wordprobs = zeros(nDirs, length(code.words));
    wordProbsT = zeros(tbins, nDirs, length(code.words));
    for d = 1:nDirs
        for w = 1:length(code.words)
            wtimes = code.code(:,:,d,:) == code.words(:, w);
            wc = sum(wtimes);
            nw = sum(wc(:) == N);
            code.wordprobs(d, w) = nw / wordcount;
            
            % For noise
            for t = 1:tbins
                wtimes = code.code(:,t,d,:) == code.words(:, w);
                wc = sum(wtimes);
                nw = sum(wc(:) == N);
                wordProbsT(t, d, w) = nw / tWC;
            end
        end
    end
    
    % Get entropies
    Pn = sum(code.wordprobs) / nDirs;
    nSsums = wordProbsT .* log2(wordProbsT);
    nSsums(isnan(nSsums)) = 0;
    noiseS = -sum(sum(sum(nSsums)));
    code.wordProbsT = wordProbsT;
    code.entropy = -sum(Pn .* log2(Pn));
    code.info = code.entropy - noiseS / (nDirs * tbins);
    
    
    code.size = length(code.words);
    code.weights = sum(code.words);
    code.sparsity = sum(code.weights) / (code.size * N);
    code.redundancy = 1 - log2(code.size) / N;
end

