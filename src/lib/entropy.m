function S = entropy(P)
%ENTROPY Calculates Shannon entropy. 0 <= P <= 1.
    S = sum(P .* log2(P));
end

