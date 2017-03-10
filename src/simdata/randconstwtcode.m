function shuffled = randconstwtcode(coding)
%RANDCONSTWTCODE Given an input neural coding (e.g. the one generated by
% the getCoding() function, generates a new set of codes by taking the
% average weight of the codewords in C, rounding this to an integer w, and
% generating constant-weight codes by randomly choosing w of the neurons to
% have an activation spike (i.e. equal 1).

    % Get properties
    tbins = coding.tbins; shuffled.tbins = tbins;
    reps = coding.reps; shuffled.reps = reps;
    length = coding.length; shuffled.length = length;
    dirs = coding.dirs;
    
    shuffled.code = zeros(length, tbins, length(dirs), reps);
    
    % Compute average weight
    w = 0;
    for t = 1:tbins
        for d = 1:length(dirs)
            for r = 1:reps
                w = w + sum(coding.code(:, t, d, r));
            end
        end
    end
    
    w = round(w / (tbins * length(dirs) * reps));
    
    shuffled.code = zeros(length, tbins, length(dirs), reps);
    
    % Get shuffled code
    for t = 1:tbins
        for d = 1:length(dirs)
            for r = 1:reps
                c = zeros(length, 1);
                c(randperm(length, w)) = 1;
                shuffled.code(:, t, d, r) = c;
            end
        end
    end
    
    % Get unique words
    unis = zeros(tbins * length(dirs) * reps, N);
    u = 1;
    
    for t = 1:tbins
        for d = 1:length(dirs)
            for r = 1:reps
                unis(u,:) = shuffled.code(:, t, d, r);
                u = u + 1;
            end
        end
    end
    
    shuffled.words = unique(unis, 'rows')';
    shuffled.size = length(shuffled.words);
    shuffled.weights = sum(shuffled.words);
    shuffled.sparsity = sum(shuffled.weights) / (shuffled.size * N);
    shuffled.redundancy = 1 - log2(shuffled.size) / N;
end
