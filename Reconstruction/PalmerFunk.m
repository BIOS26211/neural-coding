function [info] = PalmerFunk(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a function to use the... 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = dataload(n);
f = fieldnames(data);

% The data needs to be in 8ms time bins, currently it is in 2
% Since we are looking at it binarily, any point must either be 0 or 1. Any
% value 2 or higher is removed.
% It can only create words for the number of trials in the data set with
% the smallest number of trials (otherwise you lose that letter of the
% word)

s = size(data.first);

setn = zeros(s(1)/4, s(2));

datan = cell(1,n);
sizes = zeros(1,n);

for g = 1:n
    nm = f{g};
    set = data.(nm);
    s = size(set);
    for i = 1:s(1)/4
        for j = 1:s(2)
            for k = 1:s(3)
                setn(i,j,k) = set((4*i),j,k) + set(((4*i)-1),j,k) + set(((4*i)-2),j,k);
                if setn(i,j,k)>=2
                    setn(i,j,k)=1;
                end
            end
        end
    end
    datan{g} = setn;
    sizes(g) = s(3);
end

s = size(datan{1});
m = min(sizes);

words = strings(s(1),s(2),m);
for g = 1:n
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:m
                words(i,j,k) = words(i,j,k) + string(datan{g}(i,j,k));
            end
        end
    end
end

%%
% Compute the probability that there is a spike at any time in a given
% stimulus direction (P(n|t,theta), which will then be used to calculate
% the information

pnt = zeros(s(1),s(2)); %probability matrix
temp1 = zeros(1,n);
temp2 = zeros(1,n-1);
qi = zeros(s(1),s(2));
phi = zeros(s(1),s(2));
J = .16;


% Determining phi
phis = cell(1,n);
qis = cell(1,n);

for g=1:n
    for i = 1:s(1)
        for j = 1:s(2)
            for k =1:s(3)
                qi(i,j) = qi(i,j) + datan{g}(i,j,k);
            end
        end
    end
    qi = qi/s(3);
    qis{g} = qi;
    sqi = size(qi);
    for i = 1:sqi(1)
        for j = 1:sqi(2)
            if qi(i,j) == 0
                phi(i,j) = 0;
            else 
                phi(i,j) = log(qi(i,j)/(1-qi(i,j)));
            end
        end
    end
    phis{g} = phi;
end


% Calculate the pn

pfire = zeros(1,n);
for i = 1:n
    neur = datan{i};
    s = size(neur);
    temp3 = sum(neur);
    temp4 = sum(temp3);
    pfire(i) = sum(temp4);
    pfire(i) = pfire(i)/(s(1)*s(2)*s(3));
end

% Equation for probability (equation 4 in the paper)
% Sorry for the confusion, g is the i in the paper and k is the j. i and j
% here are for time and direction

s(2)= 14;

z = zeros(1,n);
for i = 1:s(1)
    for j = 1:s(2)
        for g = 1:n
            neur = datan{g};
            phi = phis{g};
            temp1(g) = phi(i,j)*sum(neur(i,j,:))/s(3);
            c=0;
            for k = 1:n
                neur2 = datan{k};
                if k~=g 
                    c=c+1;
                    temp2(c) = (J/2)*(sum(neur2(i,j,:))/s(3))*(sum(neur2(i,j,:))/s(3));
                end
            end
            z(g) = exp(-pfire(g));
            pnt(i,j) = z(g)*exp(sum(temp1)+sum(temp2));
        end
    end
end


% Calculate the maximum information 

imax = 0;
for i =1:n
    imax = imax + (pfire(i))*log2(pfire(i));
end

imax = -1*imax;

% Calculate the information 

in = zeros(s(1),s(2));

for i = 1:s(1)
    for j = 1:s(2)
        in(i,j) = -1*pnt(i,j)*log2(pnt(i,j));
    end
end

st = sum(in);

info = imax - ((1/(14*s(1)))*(sum(st)));


end
