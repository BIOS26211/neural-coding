function out = GenData(p1,p2,ndir,tbin)

% By: Emma Bsales
% University of Chicago
% February 1, 2017
% ebsales@uchicago.edu
%--------------------------------------------------------------------------
% This script generates a data set of information theory 'words' for 2
% Neurons. These data sets can either be independent or dependent of the
% direction of stimulus movement.
%--------------------------------------------------------------------------
% Inputs:   % p1 = probability of the first neuron firing
            % p2 = probability of the second neuron firing
            % ndir = number of possible direction bins
            % tbin = number of time bins
 % Outputs:     % out.type = if the data is independent ('indp') or dependent ('dept')
                % out.tBins = number of time bins
                % out.dBins = = number of possible direction bins
                % out.stimBin = direction (in 1:ndir) the stimulus is moving during the corresponding time bin
                % out.degBins = boundaries of the direction bins
                % out.pNeuron1 = probability of neuron 1 firing
                % out.pNeuron2 = probability of neuron 2 firing
                % out.fired = if the neurons fired during ea time bin (1 = yes)
                % out.words = 'word' (comb of 0s and 1s) for ea time bin
                % out.count = count for ea time bin
%--------------------------------------------------------------------------

%parameters
tBins = tbin;         % number of time bins*
dBins = ndir;          % number of direction bins; should be odd

pNeuron1 = p1;    % the probability that neuron 1 will fire in any given time bin*
pNeuron2 = p2;    % the probability that neuron 2 will fire in any given time bin*

fired = nan(2,tBins);     % create blank count vector
words = nan(tBins,2);     % create blank words vector
degBins = nan(1,dBins);       % blank vector to put the degrees of each direction in
stimBin = randsample(dBins,tBins,true);      % randomly sample from the number of direction bins for each time bin

if length(pNeuron1) ~= dBins
    pNeuron1 = repmat(pNeuron1,dBins);          % if the number of probabilities =/= the number of direction bins, fix it
    pNeuron1 = pNeuron1(1,:);
end
if length(pNeuron2) ~= dBins
    pNeuron2 = repmat(pNeuron2,dBins);
    pNeuron2 = pNeuron2(1,:);
end

stimBinSize = 180/(dBins-1);                                    % caluclate the size of the direciton bins
for i = 1:dBins
    degBins(i) = stimBinSize*i-(90+stimBinSize);              % calculate the degree of each direction
end

% start generating data
if dBins == 1
    firedN1 = [ones(1,round(tBins*(pNeuron1))) zeros(1,round(tBins*(1-pNeuron1)))];     % determine if stim indpt neuron 1 fired
    firedN2 = [ones(1,round(tBins*(pNeuron2))) zeros(1,round(tBins*(1-pNeuron2)))];     % determine if stim indpt neuron 2 fired
    
    fired(1,:) = firedN1(randperm(length(firedN1)));     % shuffle the above so its not all 1s then all 0s
    fired(2,:) = firedN2(randperm(length(firedN2)));     % shuffle the above so its not all 1s then all 0s
else
    for ii = 1:length(stimBin)
        fired(1,ii) = round(pNeuron1(stimBin(ii)));                    % determine if stim dept neuron 1 fired
        fired(2,ii) = round(pNeuron2(stimBin(ii)));                    % determine if stim dept neuron 2 fired
    end
end

for i = 1:tBins
words(i,:) = (fired(:,i))';           % put it in easier form to read as 'words'
end
count = fired(1,:)+fired(2,:);

% outputs
    if dBins == 1                        % determine if the prob is indp or dept on the direction of the stimulus
        out.type = 'indp';
    else
        out.type = 'dept';
    end
    out.tBins = tBins;
    out.dBins = dBins;
    out.stimBin = stimBin;
    out.degBins = degBins;
    out.pNeuron1 = pNeuron1;
    out.pNeuron2 = pNeuron2;
    out.fired = fired;
    out.words = words;
    out.count = count;
end
