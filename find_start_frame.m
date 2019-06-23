function [ToA] = find_start_frame(r, blockSize, CPsize)
% FIND_START_FRAME finds the index of the beginning of a frame 
% contained in signal r by computing a block based auto-correlattion 
% and applying a moving average on it (with r being a stationnary signal)
% - r is the received signal (corrupted by noise)
% - blockSize is the size of one Preamble block
% - CPsize is the Cyclic Prefix size of the Preamble blocks
 
corr = zeros(1, length(r)-2*blockSize);
for i=1:length(corr)
   corr(i)=abs(sum(r(i:i+blockSize-1).*conj(r(i+blockSize:i+2*blockSize-1))));
end
mavg = movmean(corr, [0, CPsize+1]);
[~, indx] = max(mavg);
ToA = indx - 1;
end