function hd = hammingdist(c, r)
%HAMMINGDIST Calculates the Hamming distance between vectors c and r (the
% number of positions that aren't equal).

    if (size(c) ~= size(r))
        error('Input vectors/matrices should be the same size');
    end

    hd = sum(c ~= r);
end

