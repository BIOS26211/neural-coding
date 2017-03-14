function out = GenData(n,dt,tbin,ndir,prob,stimOrder,P)

% By: Emma Bsales
% University of Chicago
% February 1, 2017. Updated March 11, 2017.
% ebsales@uchicago.edu
%--------------------------------------------------------------------------
% This script generates a data set of information theory 'words' for 'n'
% Neurons. These data sets can either be independent or dependent of the
% direction of stimulus movement.
%--------------------------------------------------------------------------
% Inputs:   % n = number of neurons in the population
            % dt = firing rate
            % prob = probability of each neuron firing for each stimulus
                % direction (should be n x ndir in size)
            % ndir = number of possible direction bins
            % tbin = number of time bins
            % P = if the data has a poisson or guassian distribution
            % stimBin = the stimulus direction bin for each time period 
                 % (must be time bins x 1 of values in 1:dBins)
               % (0 = Guass, 1 = Poiss)
 % Outputs:     % out.type = if the data is independent ('indp') or dependent ('dept')
                % out.nNeurons = number of neurons               
                % out.tBins = number of time bins
                % out.dBins = = number of possible direction bins
                % out.stimBin = direction (in 1:ndir) the stimulus is moving during the corresponding time bin
                % out.degBins = boundaries of the direction bins
                % out.probabilities = probability of neurons firing
                % out.fireRate = firing rate
                % out.distribution = the gaussian or poisson distributed
                    % random variables
                % out.fired = if the neurons fired during ea time bin (1 = yes)
                % out.words = 'word' (comb of 0s and 1s) for ea time bin
                % out.count = count for ea time bin
                % out.InputForm = put it into the form needed for the
                    % reconstruction script
                % out.InputFormData = InputForm with bin values and times
%--------------------------------------------------------------------------

% parameters
reps = 256*3;
isPoisson = P;
fireRate = dt;

tBins = tbin;         % number of time bins*
dBins = ndir;          % number of direction bins; should be odd

nNeurons = n;       % number of neurons
pNeuron = prob;        % probability of each neuron firing for each stimBin
stimOrder = stimOrder;

fired = nan(nNeurons,reps);     % create blank count vector
words = nan(reps,nNeurons);     % create blank words vector
degBins = nan(1,dBins);       % blank vector to put the degrees of each direction in
dist = nan(nNeurons,reps);        % distribution variable
stimBin = nan;
for i = 1:ndir
stimBin = [stimBin repmat(stimOrder(i),[1,reps])];
end
stimBin = stimBin(~isnan(stimBin));


if length(pNeuron(1,:)) == nNeurons
    pNeuron = pNeuron';
end
if length(pNeuron) == 1 
    pNeuron = repmat(pNeuron,[nNeurons dBins]);
elseif length(pNeuron) == dBins && length(pNeuron) ~= nNeurons
    pNeuron = repmat(pNeuron, [nNeurons 1]);
elseif length(pNeuron) ~= dBins && length(pNeuron) == nNeurons
    pNeuron = repmat(pNeuron, [1 dBins]);
end

stimBinSize = 180/(dBins-1);                                    % caluclate the size of the direciton bins
for i = 1:dBins
    degBins(i) = stimBinSize*i-(90+stimBinSize);              % calculate the degree of each direction
end

% start generating data
if isPoisson == 1
    dist = poissrnd(0.5,[nNeurons,reps]);
else
    dist = 2*normrnd(0.5,0.2,[nNeurons,reps]);
end

for jj = 1:n
for ii = 1:reps
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
sampTimes = 1:fireRate:reps;
InputFormData = zeros(length(sampTimes)+1,dBins+1,nNeurons);
for j = 1:nNeurons
    InputFormData(1,2:end,j) = degBins;
    for i = 2:length(sampTimes)+1
        st = (i-1)*fireRate;
        placement = stimBin(st-1);
        InputFormData(i,placement+1,j) = fired(j,st-1);
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
    out.stimOrder = stimOrder;
    out.degBins = degBins;
    out.probabilities = pNeuron;
    out.fireRate = fireRate;
    out.distribution = dist; % randomly generated numbers either poisson or gaussian distribution
    out.fired = fired;
    out.words = words;
    out.count = count;
    out.InputFormData = InputFormData;
    out.InputForm = InputForm;
end
