%% GDRecreation
%
% By: Emma Bsales
% University of Chicago
% March 14, 2017.
% ebsales@uchicago.edu
%--------------------------------------------------------------------------
% This script calculates the coefficient J and the peristimulus
% distributions for various generated data sets. It shows that there is not
% a difference between Gaussian and Poisson distributed data and therefore
% the information will be the same.
%--------------------------------------------------------------------------
%
%
load('probs.mat');  % load the matrix of probabilities
%
%% Part 1: Recreate J for Gaussian Data (many neurons)
% Calculate the Parameter J for Gaussian Data with all of the Neurons
% Parameters
n = length(probs); % number of neurons
dt=8; 
np = 0.15; % noise variable
nq = 0.15; % noise variable
P = 0;      % Gaussian
ndir = 14;      % number of stim directions
trials = 222;   % number of trials = largest number of trials used by Osborne
prob = probs(1:n,:);    % load only the number of probabilities we want

% Generate Gaussian Data
for j = 1:20
counts = nan(n,trials);
% Calculate counts
for i = 1:trials
stimOrder = randsample(ndir,1,true);
out = GenData(n,dt,prob,np,nq,stimOrder,P);
counts(:,i) = out.count;
end

spikes = sum(counts);

% Calculate correlation coefficients
[rho pval] = corrcoef(counts);
estJ(j) = mean(mean(rho));
end
J = max(estJ);

fprintf('Part 1:\n For Gaussian Data with all of the neurons, the estimate value of the coefficient J is %1.4f.\n \n',J)
%
%% Part 2: Generate Gaussian Data (4 neurons) with Peristimulus
% Plot Peristimulus Data for Gaussian Data with 4 Neurons
% Parameters
reps = 256*3; 
dt=8; 
np = 0.15; % noise variable
nq = 0.15; % noise variable
P = 0;      % Gaussian
ndir = 14;      % number of stim directions
trials = 222;   % number of trials = largest number of trials used by Osborne
n2 = 4; % number of neurons
prob2 = probs(1:n2,:); % load only the number of probabilities we want

% Generate Gaussian Data
fired = nan(n2,reps,trials);
for i = 1:trials
stimOrder = randsample(ndir,1,true);
out2 = GenData(n2,dt,prob2,np,nq,stimOrder,P);
fired(:,:,i) = out2.fired;
end

% Calculate frequencies
freq = periFreq(fired);
% Plot Peristimulus
corrs = peristimGraph(freq);

%% Part 3: Coefficients Gaussian vs Poisson
% Calculate the Parameter J for Gaussian and Poisson Data with all of the Neurons
% Parameters
n = length(probs); % number of neurons
dt=8; 
np = 0.15; % noise variable
nq = 0.15; % noise variable
P = 0;      % Gaussian
P2 = 1;     % Poisson
ndir = 14;      % number of stim directions
trials = 222;   % number of trials = largest number of trials used by Osborne
prob = probs(1:n,:);    % load only the number of probabilities we want

% Generate Gaussian Data
for j = 1:20
counts2 = nan(n,trials);
% Calculate counts
for i = 1:trials
stimOrder = randsample(ndir,1,true);
out = GenData(n,dt,prob,np,nq,stimOrder,P);
counts2(:,i) = out.count;
end

spikes = sum(counts2);

% Calculate correlation coefficients
[rho2 pval2] = corrcoef(counts2);
estGaussJ(j) = mean(mean(rho2));
end
J2 = max(estGaussJ);

% Generate Poisson Data
for j = 1:20
counts3 = nan(n,trials);
% Calculate counts
for i = 1:trials
stimOrder = randsample(ndir,1,true);
out = GenData(n,dt,prob,np,nq,stimOrder,P2);
counts3(:,i) = out.count;
end

spikes = sum(counts3);

% Calculate correlation coefficients
[rho3 pval3] = corrcoef(counts3);
estPoissJ(j) = mean(mean(rho3));
end
J3 = max(estPoissJ);

fprintf('Part 3:\n For Gaussian Data with all of the neurons, the estimate value of the coefficient J is %1.4f.\n For Poisson Data, the estimated value of the coefficient J is %1.4f.\n \n',J2,J3)
%
%% Part 4: Peristimulus Poisson vs Gaussian
% Plot Peristimulus Data for Gaussian and Poisson Data with 4 Neurons
% Parameters
n = length(probs); % number of neurons
dt=8; 
np = 0.15; % noise variable
nq = 0.15; % noise variable
P = 0;      % Gaussian
P2 = 1;     % Poisson
ndir = 14;      % number of stim directions
trials = 222;   % number of trials = largest number of trials used by Osborne
prob = probs(1:n,:);    % load only the number of probabilities we want
n2 = 4; 
prob2 = probs(1:n2,:);
reps = 256*3; % number of trials
fired2 = nan(n2,reps,trials);

% Calculate Gaussian Peristim
% generate gauss data
for i = 1:trials
    stimOrder = randsample(ndir,1,true);
    out2 = GenData(n2,dt,prob2,np,nq,stimOrder,P);
    fired2(:,:,i) = out2.fired;
end
% calculate frequencies
freq2 = periFreq(fired2);
% calculate peristimulus
corrs2 = peristimGraph(freq2);
% raster plot
figure;
RasterPlot(fired2(:,:,n2),'k',10,8)
  

% Calculate Poisson Peristim
% generate poiss data
for i = 1:trials
    stimOrder = randsample(ndir,1,true);
    out2 = GenData(n2,dt,prob2,np,nq,stimOrder,P2);
    fired3(:,:,i) = out2.fired;
end
% calculate frequencies
freq3 = periFreq(fired3);
% calculate peristimulus
corrs3 = peristimGraph(freq3);
% raster plot
figure;
RasterPlot(fired3(:,:,n2),'k',10,8)


 
%{
plotf = nan(2,reps);
for jj = 1:2
stimOrder = randsample(ndir,1,true);
out3 = GenData(1,dt,prob2,np,nq,stimOrder,P-jj);
plotf(jj,:) = sum(out3.fired);
time = 1:reps;
end

figure; hold on
for k = 1:2
%offset = k*max(max(plotf));
offset = 0;
plot(time,plotf(k,:)+offset,'d');
end
title('Firing Distribution across Time'); legend('Gaussian','Poisson');
xlabel('time'); ylabel('number of spikes across trials');

[r p] = corrcoef(plotf(1),plotf(2));
%}

fprintf('Part 4:\n The first two plots are of Gaussian Data with 4 neurons.\n The second two plot are of Poisson Data with 4 neurons.\n \n')

%% Part 5: Significant Testing
% Find the p-values associated with the J values found

% Check between J values calculated
propJ = mean([0.11 0.16]);
repJ = repmat(propJ,size(estJ));

[r p] = corrcoef(propJ,J);
[r2 p2] = corrcoef(propJ,J2);
[r3 p3] = corrcoef(propJ,J3);
[rJ pJ] = corrcoef(repJ,estJ);
[rG pG] = corrcoef(repJ,estGaussJ);
[rP pP] = corrcoef(repJ,estPoissJ);
[r23 p23] = corrcoef(J2,J3);
[rGP pGP] = corrcoef(estGaussJ,estPoissJ);

mins = [min(p); min(p2); min(p3); min(min(pJ)); min(min(pG)); min(min(pP)); min(p23); min(min(pGP))];
lowestPVal = min(mins);

% Use the variable pval from before
numSig = 0;
for i = 1:length(pval)
    for j = 1:length(pval)
        if pval(i,j) < 0.05
            numSig = numSig+1;
        end
    end
end
percentSig = numSig/((length(pval))^2);

numSig2 = 0;
for i = 1:length(pval2)
    for j = 1:length(pval2)
        if pval2(i,j) < 0.05
            numSig2 = numSig2+1;
        end
    end
end
percentSig2 = numSig2/((length(pval2))^2);

numSig3 = 0;
for i = 1:length(pval3)
    for j = 1:length(pval3)
        if pval3(i,j) < 0.05
            numSig3 = numSig3+1;
        end
    end
end
percentSig3 = numSig3/((length(pval3))^2);


fprintf('Part 5:\n The lowest calculated p-value between the estimated J values and the average J value from the paper is p = %1.4f.\n The percentage of p-values belowe significance for the original 3 p-value tests are p = %1.4f, p = %1.4f, and p = %1.4f.\n \n',lowestPVal,percentSig,percentSig2,percentSig3)


