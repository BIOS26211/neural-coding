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
% Roughly 5600s/trial
trials = 3;         % General neuron trials
dtrials = 10;       % Decoder trials
nneurons = 18;
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
fprintf('--------------------------------------\n\n');

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
               word = ncode.words(:, w);
               noisyw = noisegen(word, p, q);
               wML = mlDecoder(noisyw, ncode.words, p, q);
               wMAP = mapDecoder(noisyw, ncode.words, ncode.sparsity, p, q);
               
               mlCorrect = mlCorrect + sum(wML == word) / len;
               mapCorrect = mapCorrect + sum(wMAP == word) / len;
           end
        end
        mldata(n, t) = mlCorrect / (ncode.size * dtrials);
        mapdata(n, t) = mapCorrect / (ncode.size * dtrials);
    end
    
    % Get end time and calculate elapsed time
    t_elapsed = toc;
    t_total = t_total + t_elapsed;
    
    % Print info
    fprintf(' done! (elapsed time: %.3f sec)\n', t_elapsed);
    fprintf('\tInformation: %.3f\n', mean(idata(n,:)));
    fprintf('\t# redundant: %.3f\n', mean(rdata(n,:)));
    fprintf('\t%% redundant: %.3f\n', mean(prdata(n,:)));
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
plot(nn, mean(rdata(2:end), 2), 'r-');
title('Number of redundant neurons for various code sizes');
xlabel('Code size (number of neurons'); ylabel('Number of redundant neurons');

figure();
plot(nn, mean(prdata(2:end), 2), 'b-');
title('Percentage of redundant neurons for various code sizes');
xlabel('Code size (number of neurons)'); ylabel('Percentage of redundant neurons');

figure();
plot(nn, mean(idata(2:end), 2), 'g-');
title('Information of words of different neuron lengths');
xlabel('Number of neurons'); ylabel('Information (bits)');

figure();
plot(nn, mean(mldata(2:end), 2), 'c-');
title('ML decoder accuracy for various code lengths');
xlabel('Number of neurons'); ylabel('Accuracy');

figure();
plot(nn, mean(mapdata(2:end), 2), 'm-');
title('MAP decoder accuracy for various code lengths');
xlabel('Number of neurons'); ylabel('Accuracy');

%% Decoder accuracy for various noise rates
% Get a sample neuron population to work with
num_neurons = 8;
ndata = loadMTData(num_neurons);
ncode = getCoding(ndata);
len = ncode.length;
ntrials = 16;

% Define noise rates
p = 0:0.02:min(ncode.sparsity, 0.45); p(1) = 0.01;
pf = min(ncode.sparsity, 0.45) / 10;
q = p + 0.05;
qf = min(ncode.sparsity, 0.45) + 0.01;

pMLdata = zeros(length(p), ntrials);     % ML decoder data (% correct)
pMAPdata = zeros(length(p), ntrials);    % MAP decoder data (% correct)

qMLdata = zeros(length(q), ntrials);     % ML decoder data (% correct)
qMAPdata = zeros(length(q), ntrials);    % MAP decoder data (% correct)

% Testing on varying p
% Decoder trials
for i = 1:length(p)
    for d = 1:ntrials
        pMLCorrect = 0;
        pMAPCorrect = 0;
        qMLCorrect = 0;
        qMAPCorrect = 0;
        for w = 1:ncode.size
            word = ncode.words(:, w);
            % Testing on varying p
            noisyw = noisegen(w, p(i), qf);
            wML = mlDecoder(noisyw, ncode.words, p(i), qf);
            wMAP = mapDecoder(noisyw, ncode.words, ncode.sparsity, p(i), qf);

            if (sum(wML == word) == len)
                pMLCorrect = pMLCorrect + 1;
            end
            if (sum(wMAP == word) == len)
                pMAPCorrect = pMAPCorrect + 1;
            end
            
            %pMLCorrect = pMLCorrect + sum(wML == word) / len;
            %pMAPCorrect = pMAPCorrect + sum(wMAP == word) / len;

            % Testing on varying q
            noisyw = noisegen(word, pf, q(i));
            wML = mlDecoder(noisyw, ncode.words, pf, q(i));
            wMAP = mapDecoder(noisyw, ncode.words, ncode.sparsity, pf, q(i));

            if (sum(wML == word) == len)
                qMLCorrect = qMLCorrect + 1;
            end
            if (sum(wMAP == word) == len)
                qMAPCorrect = qMAPCorrect + 1;
            end
            
            %qMLCorrect = qMLCorrect + sum(wML == word) / len;
            %qMAPCorrect = qMAPCorrect + sum(wMAP == word) / len;
        end
        pMLdata(i, d) = pMLCorrect;
        pMAPdata(i, d) = pMAPCorrect;
        qMLdata(i, d) = qMLCorrect;
        qMAPdata(i, d) = qMAPCorrect;
    end
end

pMLdata = pMLdata / ncode.size;
pMAPdata = pMAPdata / ncode.size;
qMLdata = qMLdata / ncode.size;
qMAPdata = qMAPdata / ncode.size;

figure();
plot(p, mean(pMLdata, 2), 'c-');
title('ML decoder accuracy for various p rates');
xlabel('p'); ylabel('Accuracy');

figure();
plot(p, mean(pMAPdata, 2), 'm-');
title('MAP decoder accuracy for various p rates');
xlabel('p'); ylabel('Accuracy');

figure();
plot(q, mean(qMLdata, 2), 'c-');
title('ML decoder accuracy for various q rates');
xlabel('q'); ylabel('Accuracy');

figure();
plot(q, mean(qMAPdata, 2), 'm-');
title('MAP decoder accuracy for various q rates');
xlabel('q'); ylabel('Accuracy');


%% Other questions
% How do we simulate receptive fields using the data Dr. Palmer provided?
% Obviously can't do 2D RFs since the data uses 1D stimuli, but can 1D RFs
% be recreated?

% "In the case of orientation turning curves, the stimulus space is the
% interval [0, pi), and the corresponding RF code is 1D."


%%
close all;