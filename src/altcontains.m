function inds = altcontains(cells, str)
    n = length(cells);
    inds = zeros(size(cells));
    for i = 1:n
        if (~isempty(strfind(cells{i}, '.mat')))
            inds(i) = true;
        else
            inds(i) = false;
        end
    end
end