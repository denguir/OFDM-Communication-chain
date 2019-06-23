function [h] = get_ch_rsp(distrib)
%GET_CH_RSP generates channel impulse response from statistical
%distribution
% distrib is a cell array of stochastic distributions
% h is one realization of each distribution
h = zeros(size(distrib));
for i=1:length(h)
    h(i) = random(distrib{i}, 1, 1) * exp(1i*random('unif',0,2*pi));
end
end

