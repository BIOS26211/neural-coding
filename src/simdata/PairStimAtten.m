function out = PairStimAtten(n,tbin,ndir,prob,w,stimBin,P)

% By: Emma Bsales
% University of Chicago
% March 5, 2017.
% ebsales@uchicago.edu
%--------------------------------------------------------------------------
% This script generates a data set of information theory 'words' for 'n'
% Neurons. These data sets are dependent on the direction of stimulus 
% movement. The set deals with 2 stimuli per trial with the probabilities
% are independent of attention (calculated by multiplying the probability 
% (of firing for each stimulus direction).
%--------------------------------------------------------------------------
% Inputs:   % n = number of neurons in the population
            % prob = probability of each neuron firing for each stimulus
                % direction (should be n x ndir in size)
            % ndir = number of possible direction bins
            % tbin = number of time bins
            % P = if the data has a poisson or guassian distribution
                           % (0 = Guass, 1 = Poiss)
            % stimBin = the stimulus direction bin for each time period 
                 % (must be time bins x 1 of values in 1:dBins)
            % w = matrix of weights of attention to left or right stimulus
               % (2x1 matrix) weights should add up to 2 with [1 1] = equal
               % probability of choosing each side and 0.5 =more likely to
               % choose left (1 = left; 2 = right)
 % Outputs:      % out.type = if the data is independent ('indp') or dependent ('dept')
                 % out.nNeurons = number of neurons               
                 % out.tBins = number of time bins
                 % out.dBins = = number of possible direction bins
                 % out.stimBin = direction (in 1:ndir) the stimulus is moving during the corresponding time bin
                 % out.degBins = boundaries of the direction bins
                 % out.probabilities = probability of neurons firing 
                     % (calculated by multiplying the probabilities of the neuron at each of the stim bins for that time period)
                 % out.distribution = the gaussian or poisson distributed
                     % random variables
                 % out.fired = if the neurons fired during ea time bin (1 = yes)
                 % out.words = 'word' (comb of 0s and 1s) for ea time bin
                 % out.count = count for ea time bin
                 % out.InputForm = put it into the form needed for the
                     % reconstruction script
                 % out.InputFormData = InputForm with bin values and times
                 % out.attention = if they pay attention to the left or the
                     % right stimulus
                 % out.weight = weight of the left and right directions
%--------------------------------------------------------------------------

% parameters
isPoisson = P;
fireRate = 1;

tBins = tbin;         % number of time bins*
dBins = ndir;          % number of direction bins; should be odd
weights = w;
nNeurons = n;       % number of neurons
probs = prob;        % probability of each neuron firing for each stimBin
stimBin = stimBin;

fired = nan(nNeurons,tBins);     % create blank count vector
words = nan(tBins,nNeurons);     % create blank words vector
degBins = nan(1,dBins);       % blank vector to put the degrees of each direction in
%stimBin = randsample(dBins,tBins,true);      % randomly sample from the number of direction bins for each time bin
dist = nan(nNeurons,tBins);        % distribution variable

if length(probs(1,:)) == nNeurons
    probs = probs';
end
if length(probs) == 1 
    probs = repmat(probs,[nNeurons dBins]);
elseif length(probs) == dBins && length(probs) ~= nNeurons
    probs = repmat(probs, [nNeurons 1]);
elseif length(probs) ~= dBins && length(probs) == nNeurons
    probs = repmat(probs, [1 dBins]);
end

stimBinSize = 180/(dBins-1);                                    % caluclate the size of the direciton bins
for i = 1:dBins
    degBins(i) = stimBinSize*i-(90+stimBinSize);              % calculate the degree of each direction
end

% calculated probability
attention = randsample(nNeurons,tBins,true,weights);
attention = attention';

pNeuron = nan(nNeurons,tBins);
for j = 1:nNeurons
    for i = 1:tBins
        pNeuron(j,i) = probs(j,attention(i));
    end
end

% start generating data
if isPoisson == 1
    dist = poissrnd(0.5,[nNeurons,tBins]);
else
    dist = 2*normrnd(0.5,0.2,[nNeurons,tBins]);
end

for jj = 1:n
    for ii = 1:tBins
        fired(jj,ii) = round(dist(jj,ii)*pNeuron(jj,stimBin(ii)));
        if fired(jj,ii) > 1
            fired(jj,ii) = 1;
        end
    end
end
for i = 1:tBins
    words(i,:) = (fired(:,i))';           % put it in easier form to read as 'words'
end
count = sum(fired);

% fix form
InputFormData = zeros(tBins+1,dBins+1,nNeurons);
for j = 1:nNeurons
    InputFormData(1,2:end,j) = degBins;
    for i = 2:tBins+1
        placement = stimBin(i-1);
        InputFormData(i,placement+1,j) = fired(j,i-1);
        InputFormData(i,1,j) = (i-1)*fireRate;
    end
end
InputForm = InputFormData(2:end,2:end,:);

% outputs
    if dBins == 1                        % determine if the prob is indp or dept on the direction of the stimulus
        out.type = 'indp';
    else
        out.type = 'dept';
    end
    out.nNeurons = nNeurons;
    out.tBins = tBins;
    out.dBins = dBins;
    out.isPoisson = isPoisson;
    out.stimBin = stimBin;
    out.degBins = degBins;
    out.probabilities = pNeuron; % are calculated by multiplying the probabilities of the neuron at each of the stim bins for that time period
    out.distribution = dist; % randomly generated numbers either poisson or gaussian distribution
    out.fired = fired;
    out.words = words;
    out.count = count;
    out.InputFormData = InputFormData;
    out.InputForm = InputForm;
    out.attention = attention;
    out.weights = weights;

end