function pword = mapDecoder(word, codes, s, p, q)
%MAPDECORDER Given a binary word, the set of all codewords it belongs to,
% the sparsity s of the set of codewords, and the noise probabilities p (of
% a 0-bit turning into a 1) and q (of a 1-bit turning into a 0), finds the
% maximum a priori estimate of the actual codeword. p and q should naturally
% be in the range [0, 1], but no error checking for this is performed.

    map = 0;
    pword = zeros(size(word));
    for i = 1:length(codes)
        % Goes through each word in the code
        c = codes(:,i);
        wmap = sum(c .* word) * log((1 - p) * (1 - q) / (p * q)) - ...
            hammingdist(c, word) * log((1 - p) * (1 - s) / (q * s));
        if (wmap > map)
            map = wmap;
            pword = c;
        end
    end
end

