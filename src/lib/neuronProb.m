function p = neuronProb(neuron)
%NEURONPROB Given a neuron struct, calculates the trial-averaged firing
% rate probability at each time for each direction. Returns the firing
% rate probability matrix.
    
    if (~any(strcmp('data', fieldnames(neuron))))
        error('Couldnt get neuron data');
    end
    
    means = nan(36,18);
    p = means;
    for i = 1:length(neuron)
        if (max(neuron(i).data(:)) > 1)
            neuron(i).data(neuron(i).data > 1) = 1;
        end
        
        if length(mean(mean(neuron(i).data,1),3)) ~=18
            means(i,:) = [mean(mean(neuron(i).data,1),3) 0 0 0 0];
        elseif length(mean(mean(neuron(i).data,1),3)) == 18
            means(i,:) = mean(mean(neuron(i).data,1),3);
        end
    end
    
    p = 1 ./ means;
    p(isinf(p)) = 0;
end