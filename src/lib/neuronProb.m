function p = neuronProb(neuron)
%NEURONPROB Given a neuron struct, calculates the trial-averaged firing
% rate probability at each time for each direction. Returns the firing
% rate probability matrix.

    if (~any(strcmp('data', fieldnames(neuron))))
        error('Couldnt get neuron data');
    end
    
    if (max(neuron.data(:)) > 1)
        neuron.data(neuron.data > 1) = 1;
    end
    
    p = 1 ./ mean(neuron.data, 3);
    p(isinf(p)) = 0;
end

