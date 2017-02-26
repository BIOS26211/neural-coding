%% Firing vs Time Graphs
%
% Create visual representations of the firing of 2 neurons using GenData
%
% By: Emma Bsales
%% Generate the data

p1 = [1 0.95 0.85 0.80 0.75 0.65 0.55 0.50 0.45 0.35 0.25 0.15 0.05]; 
p2 = [0.05 0.15 0.25 0.35 0.45 0.50 0.55 0.65 0.75 0.80 0.85 0.95 1];

indp = GenData(0.3,0.8,1,26);
dept = GenData(p1,p2,13,26);

%% Plot the Data
% independent generated data
x = 1:indp.tBins;
y = [2*indp.fired(1,:); indp.fired(2,:)];
for j = 1:length(indp.fired)
    for jj = 1:min(size(indp.fired))
        if y(jj,j) == 0
            y(jj,j) = nan;
        end
        if y(jj,j) == 0
            y(jj,j) = nan;
        end
    end
end

figure; hold on
title('Independent of Stimulus'); xlabel('Time Fired (time bin #)'); ylabel('Neuron #');
xlim([0,indp.tBins]); ylim([0,3]);
plot(x,y,'kx','markersize',10);
set(gca,'xtick',[0:1:indp.tBins]); set(gca,'ytick',[0:1:(min(size(indp.fired)+1))]);


% dependent generated data
a = 1:dept.tBins;
b = [2*dept.fired(1,:); dept.fired(2,:)];
for j = 1:length(dept.fired)
    for jj = 1:min(size(dept.fired))
        if b(jj,j) == 0
            b(jj,j) = nan;
        end
        if b(jj,j) == 0
            b(jj,j) = nan;
        end
    end
end

figure; hold on
title('Dependent on Stimulus'); xlabel('Time Fired (time bin #)'); ylabel('Neuron #');
xlim([0,dept.tBins]); ylim([0,3]);
plot(a,b,'kx','markersize',10);
set(gca,'xtick',[0:1:dept.tBins]); set(gca,'ytick',[0:1:(min(size(dept.fired)+1))]);



