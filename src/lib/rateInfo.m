function I = rateInfo(r, dt)
%STIMINFO Given a vector r of firing rates over time, calculates the
% information provided by the rate across time. The rate can be a rate of
% any event, the occurrence of a given word, or a spike or silence from a
% particular cell.
% Parameters:
%   r: should be a single column or row vector.
%   dt: scalar number indicating the time step length in ms
% Returns:
%   I: information as a single scalar number.

    % Convert dt to ms
    dt = dt / 1000;

    % Get mean rate
    r_avg = mean(r);
    
    % Scaled rates
    r_scaled = r / r_avg;
    log_rs = log2(r_scaled); log_rs(isinf(log_rs)) = 0;
    
    I = sum(r_scaled .* log_rs .* dt) / length(r);
end

