function freq = periFreq(fired)

tl = 256*3;                 % trial length = period of interest
nt = size(fired,3);      % number of trials
nNeurons = min(size(fired));
trials = nan(nt,tl,nNeurons);
% fill trials
for i = 1:nNeurons
    for j = 1:tl
        for k = 1:nt
            trials(k,j,i) = fired(i,((k-1)*tl+j));
        end
    end
end
% sum number of spikes per time pt over all trials
nSpikes = sum(trials);

% calculate the frequency in ms
msFR = nSpikes./tl;

% change to seconds
freq = msFR*1000;

end