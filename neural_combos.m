%% Combinatoral Neural Codes from a Mathematical Coding Theory Perspective
% Work adapted from the paper by:
%   Carina Curto, Vladimir Itskov, Katherine Morrison,
%   Zachary Roth, Judy L. Walker
close all; clear all;

%% Overhead stuff -- must run before doing anything
% Add src folder to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
    addpath('src/lib');
    addpath('src/neuralcoding');
    addpath('src/simdata');
end

%% Check stats of a subset of neurons in data
% This takes a while
fprintf('Start time: %s\n', datetime('now'))
n = 6;
ndata = loadMTData(n);
ncode = getCoding(ndata);
fprintf('End time: %s\n', datetime('now'))


% Firing rate info demo
for i = 1:n
    r = neuronProb(ndata(i));
    rInfo = rateInfo(r(:, 8), 2);
    fprintf('Info from firing rate in direction 8: %.3f bits\n', rInfo);
end

%% Neuron redundancy in neural codes of various lengths
% So does this
trials = 7;
nneurons = 2:2:16; %nneurons(1) = 1;
rdata = zeros(length(nneurons), 1);
prdata = zeros(length(nneurons), 1);
idata = zeros(length(nneurons), 1);
fprintf('Start time: %s\n', datetime('now'))
fprintf('Generating...\n')
for n = 1:length(nneurons)
    for t = 1:trials
        ndata = loadMTData(nneurons(n));
        ncode = getCoding(ndata);
        r = ncode.redundancy;
        len = ncode.length;
        rdata(n) = rdata(n) + r * len;
        idata(n) = idata(n) + ncode.info;
    end
    rdata(n) = rdata(n) / trials;
    prdata(n) = rdata(n) / nneurons(n);
    idata(n) = idata(n) / trials;
    fprintf('\t%d neurons:\t%.4f redundant,\tinfo = %.4f\n', nneurons(n), rdata(n), idata(n))
end
fprintf('End time: %s\n', datetime('now'))

figure();
plot(nneurons, rdata, 'r-');
title('Number of redundant neurons for various code sizes');
xlabel('Code size (number of neurons'); ylabel('Number of redundant neurons');

figure();
plot(nneurons, prdata, 'b-');
title('Percentage of redundant neurons for various code sizes');
xlabel('Code size (number of neurons)'); ylabel('Percentage of redundant neurons');

figure();
plot(nneurons, idata, 'g-');
title('Information of words of different neuron lengths');
xlabel('Number of neurons'); ylabel('Information (bits)');

%%
close all;