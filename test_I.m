%% Test delta I for 1A x 1B for cell1 and cell2 for direction 1
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

%load data
n = loadMTData(2);
c = getCoding(n); %36x384 array of spikes and silences
data = c.code(:,:,1,:); %neuron, time bins, direction, trials
reps = c.reps; %number of trials

r1 = neuronProb(n(1));
I1 = rateInfo(r1(:,1),2); %direction 1, dt = 2ms (time bins)

r2 = neuronProb(n(1));
I2 = rateInfo(r2(:,1),2);

count1A1B = zeros(length(data),reps); %1A x 1B


for r = 1:reps
    dataA = data(1,:,r); %row in data cell corresponding to cell A spike train
    dataB = data(2,:,r); %row in data cell corresponding to cell B spike train
    spikeA = find(dataA == 1); %location of 1's in A
    spikeB = find(dataB == 1); %location of 1's in B
    for j = 1:length(spikeA)
        valA = spikeA(j); %time position of spike A
        for k = 1:length(spikeB)
            valB = spikeB(k); %time position of nonspike B
            diff = abs(valA - valB); %time(bins) between each spike
            if diff <= 5 %within 10 ms of each other
                count1A1B(spikeA(j)) = count1A1B(spikeA(j)) + 1;
            end
        end
    end
end
count1A1B = count1A1B ./ reps;

I12 = rateInfo(count1A1B,2);

deltaI = I12 - I1 - I2


