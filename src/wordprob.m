function prob = wordprob(W, P)
%WORDPROB Given a word W (vertical vector of 0s and 1s) and firing rates P
%  (vertical vector of the same length as W), calculates the probability of
%  the word occuring at a given point in time.
    if (length(W) ~= length(P))
        error('Vectors W and P should both be of the same length');
    end
    prob = -sum(W .* P .* log2(W .* P));
end

