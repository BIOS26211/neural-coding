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
ndata = loadMTData(2);
ncode = getCoding(ndata);

%% Neuron redundancy in neural codes of various lengths
% So does this
trials = 7;
nneurons = 0:3:36; nneurons(1) = 1;
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
    fprintf('\t%d neurons: %d redundant, info = %f\n', nneurons(n), rdata(n), idata(n))
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