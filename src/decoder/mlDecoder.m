function pword = mlDecoder(word, codes, p, q)
%MLDECODER Given a binary word, the set of all codewords it belongs to, and
% the noise probabilities p (of a 0-bit turning into a 1) and q (of a 1-bit
% turning into a 0), finds the maximum likelihood estimate of the actual
% codeword. p and q should naturally be in the range [0, 1], but no error
% checking for this is performed.

    ml = nan;
    pword = zeros(size(word));
    if (p ~= q)
        for i = 1:length(codes)
            % Goes through each word in the code
            c = codes(:,i);
            wml = sum(c .* word) * log((1 - p) * (1 - q) / (p * q)) - ...
                sum(c) * log((1 - p) / q);
            if (wml > ml || isnan(ml))
                ml = wml;
                pword = c;
            end
        end
    else
        % If p == q, minimize the Hamming distance
        for i = 1:length(code)
            c = codes(:, i);
            mld = hammingdist(c, word);
            if (mld < ml || isnan(ml))
                ml = mld;
                pword = c;
            end
        end
    end
end

