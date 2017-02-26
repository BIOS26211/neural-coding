%% Recreation - Osborne/Palmer 2008
%
% Garrett Healy
% 
% 25FEB2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to use information theory to recreate the data
% modification completed by Stephanie Palmer. These data are from MT neurons
% in macaques (n=3) anesthetized with sufentanil. 
% 
% DATA -------------------------------------------------------------------
% 
% There are 13 directions + 1 null direction. 
% 
% Recordings made for 256ms, put into 2ms time bins for 3 monkeys. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 

% Creates a raster of one set of data, to visualize it. 

cell1 = load('cell_11l04.mat');

data1 = cell1.data;
dirs = cell1.dirs;

s=size(data1);

figure 
hold on 
xlabel('Time(ms)');ylabel('Direction');
title('Neural Spikes by Stimulus Direction');
for i = 1:s(1)
    for j= 1:s(2)-1
        for k = 1:s(3)
            if data1(i,j,k) == 1
                spkx = [i*2, i*2];
                spky = [dirs(j) - 5,dirs(j) + 5];
                line(spkx,spky,'color','k','LineWidth',1);
            else 
            end
        end
    end
end
hold off
%%
% Load 4 more cells in order to make a 'word'

cell2 = load('cell_11l07.mat');
data2 = cell2.data;
cell3 = load('cell_11l10.mat');
data3 = cell3.data;
cell4 = load('cell_11l14.mat');
data4 = cell4.data;
cell5 = load('cell_11l15.mat');
data5 = cell5.data;

% The data needs to be in 8ms time bins, currently it is in 2
% Since we are looking at it binarily, any point must either be 0 or 1. Any
% value 2 or higher is removed. 

data1n = zeros(s(1)/4, s(2));
data2n = zeros(s(1)/4, s(2));
data3n = zeros(s(1)/4, s(2));
data4n = zeros(s(1)/4, s(2));
data5n = zeros(s(1)/4, s(2));

for i = 1:s(1)/4
    for j = 1:s(2)
        data1n(i,j) = data1((4*i),j) + data1(((4*i)-1),j) + data1(((4*i)-2),j);
        if data1n(i,j)>=2
            data1n(i,j)=1;
        end
        data2n(i,j) = data2((4*i),j) + data2(((4*i)-1),j) + data2(((4*i)-2),j);
        if data2n(i,j)>=2
            data2n(i,j)=1;
        end
        data3n(i,j) = data3((4*i),j) + data3(((4*i)-1),j) + data3(((4*i)-2),j);
        if data3n(i,j)>=2
            data3n(i,j)=1;
        end
        data4n(i,j) = data4((4*i),j) + data4(((4*i)-1),j) + data4(((4*i)-2),j);
        if data4n(i,j)>=2
            data4n(i,j)=1;
        end
        data5n(i,j) = data5((4*i),j) + data5(((4*i)-1),j) + data5(((4*i)-2),j);
        if data5n(i,j)>=2
            data5n(i,j)=1;
        end
    end
end

s = size(data1n);

words = strings(s(1),s(2));
for i = 1:s(1)
    for j = 1:s(2)
        words(i,j) = string(data1n(i,j,1))+string(data2n(i,j,1))+string(data3n(i,j,1))+ ...
            string(data4n(i,j,1))+string(data5n(i,j,1));
    end
end
%%
% Compute the probability that there is a spike at any time in a given
% stimulus direction 

pnt = zeros(s(1),1);
temp1 = zeros(1,s(2)-1);
temp2 = zeros(1,s(2)-1);
J = .13;
dJ = .001*J;

for i = 1:s(2)
    for k = 1:s(2)
        
        for a = 1:s(2)
            for b = 1:s(2)-1
                if a~=b
                    
                end
            end
        end
    end
    sum(
end
phi = dJ/2;

for i = 1:s(1)
    for j = 1:s(2)
        temp1(j) = phi*data1n(i,j);
        for k = 1:s(2)-1
            if j ~= k
                temp2(k) = (J/2)*data1n(i,j)*data1n(i,k);
            end
        end
    end
    pnt(i) = exp(sum(temp1)+sum(temp2));
end

m = max(pnt);
pnt = pnt/m;
%%
% Compute the mutual information

imax =-1*(5/2)*log2(1/2);
iwords = zeros(s(1),1);

for i = 1:s(1)
    for j = 1:s(2)
        itemp = -1*pnt(i)*log2(pnt(i));
        iwords(i) = imax - ((1/(13*s(2)))*itemp);
    end
end


%%
vals = zeros(s(1),s(2)); 
trials = zeros(s(1),s(3));
spikes = 0:23;
d = cumsum(data1);
tots = zeros(1,s(2));
newinf = zeros(s(1),1);
stdvals = zeros(s(1),1);

for q = 1:s(1)
    counts = zeros(length(spikes),s(2));
    unc = zeros(1,s(2));
    for i = 1:s(2)
        for j = 1:s(3)
           x = d(q,i,j); 
           counts(x+1,i) = counts(x+1,i) + 1;
        end
    end
    for i = 1:s(2)
        for j = 1:length(spikes)
            tot2 = sum(counts(j,:));
            if tot2 ~=0
                prob = counts(j,i)/ 184;
                mprob = mean(counts(j,:))/184;
                if prob~=0
                    unc(i) = unc(i) - prob*log2(prob/mprob);
                end
            end
        end
    end
    newinf(q) = - mean(unc);
    stdvals(q) = std(unc);
end

figure
hold on
%plot(time*1000,newinf)
errorbar(time*1000,newinf,stdvals)
xlabel('Time(ms)');ylabel('Information(bits)');
title('Information Contained in Spike Counts with Regard to Stimulus Direction');
hold off





