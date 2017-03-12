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
% Note to self: Loading all neurons into memory takes > 1 hr on Surface.
fprintf('Start time: %s\n', datetime('now'))
n = 16;
ndata = loadMTData(n);
ncode = getCoding(ndata);
scode = shuffledcode(ncode);
cwtcode = constwtcode(ncode);
fprintf('End time: %s\n', datetime('now'))

%% Neuron redundancy in neural codes of various lengths
% So does this
trials = 1;
nneurons = 2:2:16;
rdata = zeros(length(nneurons), 1);
prdata = zeros(length(nneurons), 1);
idata = zeros(length(nneurons), 1);
fprintf('Generating %d data points....\n|', length(nneurons))
tic;
for n = 1:length(nneurons)
    for t = 1:trials
        ndata = loadMTData(nneurons(n));
        ncode = getCoding(ndata);
        r = ncode.redundancy;
        len = ncode.length;
        rdata(n) = rdata(n) + r * len;
        idata(n) = idata(n) + ncode.info;
        fprintf('~')
    end
    rdata(n) = rdata(n) / trials;
    prdata(n) = rdata(n) / nneurons(n);
    idata(n) = idata(n) / trials;
    %fprintf('\t%d neurons:\t%.4f redundant,\tinfo = %.4f\n', nneurons(n), rdata(n), idata(n))
    fprintf('=')
end
t_elapsed = toc;
fprintf('Elapsed time: %s s\n', t_elapsed)

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

% Try: loading up n neurons and defining that as an RF of size n. Need to
% define a binary response map (p. 1895 [5/35]) that maps a stimulus to a
% code (set of codewords).

%%
close all;