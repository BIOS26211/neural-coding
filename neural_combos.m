%% Combinatoral Neural Codes from a Mathematical Coding Theory Perspective
% Work adapted from the paper by:
%   Carina Curto, Vladimir Itskov, Katherine Morrison,
%   Zachary Roth, Judy L. Walker
close all; clear all;

% Add src folder to path
if (isempty(strfind(pwd(), strcat(filesep, 'src'))))
    addpath('src');
    addpath('src/neuralcoding');
end

trials = 10;
nneurons = 0:3:30; nneurons(1) = 1;
rdata = zeros(length(nneurons), 1);
prdata = zeros(length(nneurons), 1);
fprintf('Generating...')
for n = 1:length(nneurons)
    for t = 1:trials
        ndata = loadMTData(nneurons(n));
        ncode = getCoding(ndata);
        r = ncode.redundancy;
        len = ncode.length;
        rdata(n) = rdata(n) + r * len;
    end
    rdata(n) = rdata(n) / trials;
    prdata(n) = rdata(n) / nneurons(n);
    fprintf('%d neurons: %d redundant\n', nneurons(n), rdata(n))
end

figure();
plot(nneurons, rdata, 'r-');
title('Number of redundant neurons for various code sizes');
xlabel('Code size (number of neurons'); ylabel('Number of redundant neurons');

figure();
plot(nneurons, prdata, 'b-');
title('Percentage of redundant neurons for various code sizes');
xlabel('Code size (number of neurons)'); ylabel('Percentage of redundant neurons');
