function p = neuronProb(neuron)
%NEURONPROB Given a neuron struct, calculates the trial-averaged firing
% rate probability at each time for each direction. Returns the firing
% rate probability matrix.

    if (~any(strcmp('data', fieldnames(neuron))))
        error('Couldnt get neuron data');
    end
    
    for i = 1:length(neuron)
        if (max(neuron(i).data(:)) > 1)
            neuron(i).data(neuron(i).data > 1) = 1;
        end
    
        p(:,:,i) = 1 ./ mean(neuron(i).data, 3);
    end
    p(isinf(p)) = 0;
        
    neuron.firerate = p;
    
end