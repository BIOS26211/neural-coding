%% Combinatoral Neural Codes from a Mathematical Coding Theory Perspective
% Work adapted from the paper by:
%   Carina Curto, Vladimir Itskov, Katherine Morrison,
%   Zachary Roth, Judy L. Walker
close all; clear all;

%% Overhead stuff -- must run before doing anything
% Add src folder to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
    addpath('src/decoder');
    addpath('src/lib');
    addpath('src/neuralcoding');
    addpath('src/simdata');
end

%% Check stats of a subset of neurons in data
% Note to self: Loading all neurons into memory takes > 1 hr on Surface.
tic;
n = 8;
ndata = loadMTData(n);
ncode = getCoding(ndata);
scode = shuffledcode(ncode);
cwtcode = constwtcode(ncode);
t_elapsed = toc;
fprintf('Elapsed time: %.3f s\n', t_elapsed)

%% Neuron redundancy in neural codes of various lengths
% So does this
trials = 1;
nneurons = 2:2:18;
rdata = zeros(length(nneurons), 1);
prdata = zeros(length(nneurons), 1);
idata = zeros(length(nneurons), 1);
fprintf('Generating %d data points....\n', nneurons(end))
t_total = 0;
for n = 1:length(nneurons)
    tic;
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
    t_elapsed = toc;
    t_total = t_total + t_elapsed;
    fprintf('%d neurons:\t%.4f redundant,\tinfo = %.4f (elapsed time: %.3f s)\n', nneurons(n), rdata(n), idata(n), t_elapsed)
end
fprintf('Total elapsed time: %.3f s\n', t_total);

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

%% Other questions
% How do we simulate receptive fields using the data Dr. Palmer provided?
% Obviously can't do 2D RFs since the data uses 1D stimuli, but can 1D RFs
% be recreated?

% "In the case of orientation turning curves, the stimulus space is the
% interval [0, pi), and the corresponding RF code is 1D."


%%
close all;