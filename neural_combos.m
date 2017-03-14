%% Combinatoral Neural Codes from a Mathematical Coding Theory Perspective
% Work adapted from the paper by:
%   Carina Curto, Vladimir Itskov, Katherine Morrison,
%   Zachary Roth, Judy L. Walker
close all; clear all;

%% Overhead stuff -- must run before doing anything
% Add src folder and children to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
    addpath('src/decoder');
    addpath('src/lib');
    addpath('src/neuralcoding');
    addpath('src/simdata');
end

%% Neuron redundancy in neural codes of various lengths
% Meta
trials = 1;         % General neuron trials
dtrials = 10;       % Decoder trials
nneurons = 16;
nstep = 1;

% Parameters
p = 0.06; q = 0.1;

% Run over:
nn = 2:nstep:nneurons;

% Initialize data matrices
rdata = zeros(length(nn), trials);      % redundancy
prdata = zeros(length(nn), trials);     % percent redundant
idata = zeros(length(nn), trials);      % information
szdata = zeros(length(nn), trials);     % code set size
spdata = zeros(length(nn), trials);     % sparsity
mldata = zeros(length(nn), trials);     % ML decoder data (% correct)
mapdata = zeros(length(nn), trials);    % MAP decoder data (% correct)

fprintf('Generating %d data points....\n', length(nn))

% Run and time
t_total = 0;
for n = 2:nstep:nneurons
    fprintf('Working on %d-neuron words... ', n)
    
    % Run 
    tic;
    for t = 1:trials
        ndata = loadMTData(n);
        ncode = getCoding(ndata);
        r = ncode.redundancy;
        len = ncode.length;
        rdata(n, t) = r * len;
        prdata(n, t) = (r * len) / n;
        idata(n, t) = ncode.info;
        szdata(n, t) = ncode.size;
        spdata(n, t) = ncode.sparsity;
        
        % Decoder trials
        mlCorrect = 0;
        mapCorrect = 0;
        for d = 1:dtrials
           for w = 1:ncode.size
               noisyw = noisegen(w, p, q);
               wML = mlDecoder(noisyw, ncode.words, p, q);
               wMAP = mapDecoder(noisyw, ncode.words, ncode.sparsity, p, q);
               
               mlCorrect = mlCorrect + sum(wML == w) / len;
               mapCorrect = mapCorrect + sum(wMAP == w) / len;
           end
        end
        mldata(n, t) = mlCorrect / (w * dtrials);
        mapdata(n, t) = mapCorrect / (w * dtrials);
    end
    
    % Get end time and calculate elapsed time
    t_elapsed = toc;
    t_total = t_total + t_elapsed;
    
    % Print info
    fprintf(' Done! (Elapsed time: %.3f sec\n', t_elapsed);
    fprintf('\tInformation: %.3f\n', mean(idata(n,:)));
    fprintf('\t# redundant: %.3f\n', mean(rdata(n,:)));
    fprintf('\t% redundant: %.3f\n', mean(prdata(n,:)));
    fprintf('\tML decoder accuracy: %.3f\n', mean(mldata(n,:)));
    fprintf('\tMAP decoder accuracy: %.3f\n\n', mean(mapdata(n,:)));
end
fprintf('Total elapsed time: %.3f sec\n', t_total);

%% Figure plotting
% For reference
% rdata = zeros(length(nn), trials);      % redundancy
% prdata = zeros(length(nn), trials);     % percent redundant
% idata = zeros(length(nn), trials);      % information
% szdata = zeros(length(nn), trials);     % code set size
% spdata = zeros(length(nn), trials);     % sparsity
% mldata = zeros(length(nn), trials);     % ML decoder data (% correct)
% mapdata = zeros(length(nn), trials);    % MAP decoder data (% correct)

figure();
plot(nn, rdata, 'r-');
title('Number of redundant neurons for various code sizes');
xlabel('Code size (number of neurons'); ylabel('Number of redundant neurons');

figure();
plot(nn, prdata, 'b-');
title('Percentage of redundant neurons for various code sizes');
xlabel('Code size (number of neurons)'); ylabel('Percentage of redundant neurons');

figure();
plot(nn, idata, 'g-');
title('Information of words of different neuron lengths');
xlabel('Number of neurons'); ylabel('Information (bits)');

figure();
plot(nn, mldata, 'c-');
title('ML decoder accuracy for various code lengths');
xlabel('Number of neurons'); ylabel('Accuracy');

figure();
plot(nn, mapdata, 'm-');
title('MAP decoder accuracy for various code lengths');
xlabel('Number of neurons'); ylabel('Accuracy');

%% Other questions
% How do we simulate receptive fields using the data Dr. Palmer provided?
% Obviously can't do 2D RFs since the data uses 1D stimuli, but can 1D RFs
% be recreated?

% "In the case of orientation turning curves, the stimulus space is the
% interval [0, pi), and the corresponding RF code is 1D."


%%
close all;