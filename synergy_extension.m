clear all;
close all; 

% Add src folder to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
    addpath('src/neuralcoding');
    addpath('src/lib');
    addpath('Reconstruction');
    addpath('MT_data');
end

%% load data
n = loadMTData(20);
c = getCoding(n); %36x384 array of spikes and silences
data = c.code(:,:,8,1);
reps = c.reps;

%% 1A x 1B
%matrix of possible combinations for pairs of A and B
b2 = nchoosek(1:length(n),2); %630 possible combinations for pairs of cells

%initialize delta I array
count1A1B = zeros(length(data),length(b2)); %1A x 1B

%1A x 1B
for r = 1:reps
    for i = 1:length(b2)
        cellA = b2(i,1); %cell A index
        cellB = b2(i,2); %cell B index
        dataA = data(cellA,:); %row in data cell corresponding to cell A spike train
        dataB = data(cellB,:); %row in data cell corresponding to cell B spike train  
        spikeA = find(dataA == 1); %location of 1's in A
        spikeB = find(dataB == 1); %location of 1's in B
        for j = 1:length(spikeA)
            valA = spikeA(j); %time position of spike A
            for k = 1:length(spikeB)
                valB = spikeB(k); %time position of nonspike B
                diff = abs(valA - valB); %time(bins) between each spike
                if diff <= 5 %within 10 ms of each other
                    count1A1B(spikeA(j),i) = count1A1B(spikeA(j), i) + 1;
                end
            end
        end
    end
end
count1A1B = count1A1B ./ reps;

nI = zeros(size(data, 1), 1);
for i = 1:length(n)
    nI(i) = rateInfo(neuronProb(n(i)), 2);
end

IAB = zeros(length(b2),1);
deltaI = zeros(length(b2), 1);
for i = 1:length(b2)
    IAB(i) = rateInfo(count1A1B(:,i),2);
    if isnan(IAB(i))
        IAB(i) = 0;
    end
    deltaI(i) = IAB(i) - nI(b2(i, 1)) - nI(b2(i, 2));
end   

edges = min(deltaI):5e-04:max(deltaI);

figure;
hold on;
histogram(deltaI,edges);
xlabel('synergy (bits)');
ylabel('Probability Density');
title('Probability Density for 1A x 1B');
hold off;
    
%% 1A x 10
%initialize delta I array
count1A0B = zeros(1,length(b2)); %1A x 0B

%1A x 0B
for i = 1:length(b2)
    cellA = b2(i,1); %cell A index
    cellB = b2(i,2); %cell B index
    dataA = data(cellA,:); %row in data cell corresponding to cell A spike train
    dataB = data(cellB,:); %row in data cell corresponding to cell B spike train  
    spikeA = find(dataA == 1); %location of 1's in A
    nonspikeB = find(dataB == 0); %location of 0's in B
    for j = 1:length(spikeA)
        valA = spikeA(j); %time position of spike A
        for k = 1:length(nonspikeB)
            valB = nonspikeB(k); %time position of nonspike B
            diff = abs(valA - valB); %time(bins) between each spike
            if diff <= 25 %within 50 ms of each other
                count1A0B(1,i) = count1A0B(1,i)+1;
            end
        end
    end
end

count1A0B = count1A0B ./ reps;

nI = zeros(size(data, 1), 1);
n0 = zeros(size(data, 1), 1);
for i = 1:length(n)
    nI(i) = rateInfo(neuronProb(n(i)), 2);
    n0(i) = rateInfo(neuronProb(n(i))./(neuronProb(n(i))-1),2);
end

IAB = zeros(length(b2),1);
deltaI = zeros(length(b2), 1);
for i = 1:length(b2)
    IAB(i) = rateInfo(count1A1B(:,i),2);
    if isnan(IAB(i))
        IAB(i) = 0;
    end
    deltaI(i) = IAB(i) - nI(b2(i, 1)) - n0(b2(i, 2));
end

edges = min(deltaI):5e-04:max(deltaI);

figure;
hold on;
histogram(deltaI,edges);
xlabel('synergy (bits)');
ylabel('Probability Density');
title('Probability Density for 1A x 0B');
hold off;

%% 1A x 0B x 0C

%matrix of possible combinations for triplets of A, B, and C
b3 = nchoosek(1:length(n),3); %7140 possible combinations for triplets of cells

%initialize delta I array
count1A0B0C = zeros(1,length(b3)); %1A x 0B x 0C
count0B0C = zeros(1,length(b3));

%1A x 0B x 0C
for r = 1:reps
for i = 1:length(b3)
    cellA = b3(i,1); %cell A index
    cellB = b3(i,2); %cell B index
    cellC = b3(i,3); %cell C index
    dataA = data(cellA,:); %row in data cell corresponding to cell A spike train
    dataB = data(cellB,:); %row in data cell corresponding to cell B spike train  
    dataC = data(cellC,:); %row in data cell corresponding to cell C spike train 
    spikeA = find(dataA == 1); %location of 1's in A
    nonspikeB = find(dataB == 0); %location of 0's in B
    nonspikeC = find(dataC == 0); %location of 0's in C
    for j = 1:length(spikeA)
        valA = spikeA(j); %time position of spike A
        for k = 1:length(nonspikeB)
            valB = nonspikeB(k); %time position of nonspike B
            diffAB = abs(valA - valB); %timebins between each spike
            if diffAB <= 25 %within 50 ms of each other
                for l = 1:length(nonspikeC)
                    valC = nonspikeC(l); %time position of nonspike C
                    diffBC = abs(valB - valC);
                    if diffBC <= 10 %within 20 ms of each other
                        count1A0B0C(1,i) = count1A0B0C(1,i)+1;
                    end
                end
            end
        end
    end
    for m = 1:length(nonspikeB)
        valB0 = nonspikeB(m);
        for u = 1:length(nonspikeC)
            valC0 = nonspikeC(u);
            diffB0C0 = abs(valB0 - valC0);
            if diffB0C0 <= 10
                count0B0C(1,i) = count0B0C(1,i)+1;
            end
        end
    end        
end
end

count1A0B0C = count1A0B0C ./ reps;
count0B0C = count0B0C ./ reps;

nI = zeros(size(data, 1), 1);
for i = 1:length(n)
    nI(i) = rateInfo(neuronProb(n(i)), 2);
end

IABC = zeros(length(b3),1);
I0BC = zeros(length(b3),1);
deltaI = zeros(length(b3), 1);
for i = 1:length(b3)
    IABC(i) = rateInfo(count1A0B0C(:,i),2);
    if isnan(IABC(i))
        IABC(i) = 0;
    end
    I0BC(i)= rateInfo(count0B0C(:,i),2);
    if isnan(I0BC(i))
        I0BC(i) = 0;
    end
    deltaI(i) = IABC(i) - nI(b3(i, 1)) - I0BC(i);
end  

edges = min(deltaI):5e-05:max(deltaI);

figure;
hold on;
histogram(deltaI,edges);
xlabel('synergy (bits)');
ylabel('Probability Density');
title('Probability Density for 1A x 0B x 0B');
hold off;
