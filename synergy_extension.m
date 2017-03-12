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
%%

%load data
n = loadMTData(36);
c = getCoding(n); %36x384 array of spikes and silences
reps = c.reps;
data = c.code(:,:,1,:);

%% 1A x 1B
%matrix of possible combinations for pairs of A and B
b2 = nchoosek(1:36,2); %630 possible combinations for pairs of cells

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

%% 1A x 0B x 0C

%matrix of possible combinations for triplets of A, B, and C
b3 = nchoosek(1:36,3); %7140 possible combinations for triplets of cells

%initialize delta I array
count1A0B0C = zeros(1,length(b3)); %1A x 0B x 0C

%1A x 0B x 0C
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
end


%% test data
data25 = data(25,:);
data26 = data(26,:);
spike25 = find(data25 == 1);
spike26 = find(data26 == 1);
for i = 1:length(spike25)
    k = 0;
    val25 = spike25(i);
    for j = 1:length(spike26)
        val26 = spike26(j);
        diff = abs(val25 - val26);
        if diff <= 5
            k = k+1;
        end
    end
end
        
