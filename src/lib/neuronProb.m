function p = neuronProb(neuron)
%NEURONPROB Given a neuron struct, calculates the trial-averaged firing
% rate probability at each time 1 for each direction. Appends the info to
% the struct as well as returns the firing rate probability matrix.

    if (~any(strcmp('data', fieldnames(neuron))))
        error('Couldnt get neuron data');
    end
    
    if (max(neuron.data(:)) > 1)
        neuron.data(neuron.data > 1) = 1;
    end
    
    p = mean(neuron.data, 3);
    neuron.firerate = p;

end

