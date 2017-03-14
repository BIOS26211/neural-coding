function vec = noisegen(vec, p, q)
%NOISEGEN Takes a binary vector and, for each element, adds noise. Each 0
% is turned into a 1 with probability p, and each 1 is turned into a 0 with
% probability q. If either p or q are greater than 1, they're taken as 1.
% Typically, p <= q < 1/2, and p < s (if known) < 1/2

    % Deal with 0's
    iZeros = find(vec == 0);
    pZeros = rand(size(iZeros));
    
    % Deal with 1's
    iOnes = find(vec == 1);
    pOnes = rand(size(iOnes));

    % Do noise
    vec(iZeros & pZeros < p) = 1;
    vec(iOnes & pOnes < q) = 0;
end

