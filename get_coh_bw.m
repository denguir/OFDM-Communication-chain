function [coh_bw, delay_spread] = get_coh_bw(tau, PDP)
%GET_COH_BW computes the coherence bandwidth from a given PDP
dt = tau(2) - tau(1);
total_power = sum(PDP) * dt;
mean_tau = sum(tau .* PDP) * dt / total_power;
delay_spread = sqrt((1/total_power) * (sum((tau.^2) .* PDP * dt)) - mean_tau^2 );
coh_bw = 1/(2*pi*delay_spread);
end

