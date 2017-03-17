function corrs = peristimGraph(freq)

nNeurons = length(freq(1,1,:));
offset = round(max(max(freq(1,:,:))));
reps = length(freq);
time = 1:reps;

grD = nan(size(freq));
for i = 1:nNeurons
    grD(1,:,i) = freq(1,:,i)+(offset*(i-1));
end

figure; hold on
for j = 1:nNeurons
plot(time,grD(1,:,j))
end
title('Peristimulus');
xlabel('time (ms)'); ylabel('frequency (Hz)');

corrs = grD;

end