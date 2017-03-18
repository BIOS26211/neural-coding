function [probmatrix] = probFind(n)

x = dataload(n);
f = fieldnames(x); 
l = length(f);
probmatrix = zeros(l,14);

for i = 1:l
    n = f{i};
    neur = x.(n);
    s = size(neur);
    for h = 1:14
        for j = 1:s(1)
            for k = 1:s(3)
                probmatrix(i,h) = probmatrix(i,h) + neur(j,h,k);
            end
        end
    end
    probmatrix(i,:) = probmatrix(i,:)./(s(1)*s(3));
end
end