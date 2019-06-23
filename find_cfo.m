function [CFO] = find_cfo(r, blockSize, CPsize, ts)
%FIND_CFO Summary of this function goes here
%   Detailed explanation goes here

corr = 0;
preamble = r(CPsize+1:CPsize+2*blockSize);
for n=1:length(preamble)-blockSize
    corr = corr + preamble(n) * conj(preamble(n + blockSize));
end
CFO = angle(corr)/(2*pi*blockSize*ts);