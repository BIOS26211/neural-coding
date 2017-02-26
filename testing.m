% Testing some of the functions out
close all; clear all;
% Add src folder to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
end

% Load up all of the recorded neuron data
neuronData = loadMTData();
numNeurons = size(neuronData, 2);
raster = zeros(size(neuronData(1).data, 1), numNeurons);

dt = 0.002;
for i = 1:numNeurons
    numT = numel(neuronData(i).data(1, 1, :));
    spike_counts = sum(neuronData(i).data, 3);
    
    % Firing rate given direction and time
    neuronData(i).stimdepfirerate = spike_counts ./ (numT * dt);
    
    % Firing rate given time
    neuronData(i).firerate = sum(neuronData(i).stimdepfirerate, 2);
    raster(:,i) = squeeze(neuronData(i).data(:,8,1));
end
raster(raster(:) ~= 0) = 1;  % fix the binarization
figure; imagesc(raster'); colormap(gray);