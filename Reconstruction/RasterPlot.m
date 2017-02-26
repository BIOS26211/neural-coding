function RasterPlot(Spikes,color,fs,dt) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function plots a rasterplot
% 
% INPUTS:
% Spikes - can by spike time cell or binary matrix with dim ncells x tp
%            We assume all spike times are in ms!
% color  - plot color
% fs     - fontsize
% dt     - time step (ms) (only used for binary spike matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We don't need dt if Spikes is cell structure of spike times
if nargin == 3
    dt = [];
end

hold all

% if the spikes are in a struct, this will access that (work in progress)
% if isstruct(Spikes)
%     fuckthis = fieldnames(Spikes);
%     name = string(fuckthis{1,1});
%     data = getfield(Spikes, name);
%     L = structfun(@(field) length(field),Spikes);
%     for ii = 1:L
%         trial = Spikes.name{1,ii};
%         for jj = 1:length(trial)
%             spkx = [trial(jj), trial(jj)];
%             spky = [ii - 0.4,ii + 0.4];
%             line(spkx,spky,'color',color,'LineWidth',1);
%         end
%     end

% if spike times, then use spike times as x-values
if iscell(Spikes)
    l = length(Spikes);
    for i = 1:l
        trial = Spikes{:,i};
        for q = 1:length(trial)
            spkx = [trial(q), trial(q)];
            spky = [i - 0.4,i + 0.4];
            line(spkx,spky,'color',color,'LineWidth',1);
        end
    end
% if Spikes is not a cell structure, it is a binary matrix. Use non-zero 
% values multiplied by dt as x-values
else
    for i = 1:size(Spikes,1)
        J = find(Spikes(i,:));
        for q = 1:length(J)
            spkx = [J(q),J(q)] .* dt;
            spky = [i - 0.4,i + 0.4];
            line(spkx,spky,'color',color,'LineWidth',1);
        end
    end
end

set(gca,'fontsize',fs)
xlabel('Time (s)','fontsize',fs);
ylabel('neuron','fontsize',fs);
title('Raster');

hold off

end