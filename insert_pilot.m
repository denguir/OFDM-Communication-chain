function [bits] = insert_pilot(bits, pilotPos, Nbps, Nchannels)
%INSERT_PILOT inserts pilot bits to the pilotPos indeces

Nbits = length(bits);
pilot_bits = randi([0,1], Nbps, 1);

numBlocks = Nbits/(Nbps*Nchannels);
for i=0:numBlocks-1
    for k=pilotPos
        n1 = (Nchannels*i+k-1)*Nbps+1;
        n2 = (Nchannels*i+k)*Nbps;
        bits(n1:n2) = pilot_bits;
    end
end

