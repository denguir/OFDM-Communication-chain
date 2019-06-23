function [ data_out ] = demapping(data_in, Nbps)
%
% INPUTS:
%   - input symbol streams for each user [NbSymbols x Number of subcarriers]
%		[vector with indexing (NbOfBlocks, B)]
%   - scaling_soft gives local relative noise/interference power [Ntx x B]
%
% OUTPUTS:
%   - Demapped bit streams for each user [(NbSymbols*Nbps) x K]
%		[vector with indexing (NbOfBlocks, Ntx, B, Nbps)]


NbSymbols = size(data_in, 1);
K = size(data_in, 2);

Kmod = sqrt([1,2,NaN,10,NaN,42,NaN,170]);   % Constellation scaling factors
data_in = data_in * Kmod(Nbps);     % Undo Constellation scaling

data_out = zeros(NbSymbols*Nbps,K);

if Nbps == 1
    % BPSK demodulation
    data_out = real(data_in);
else
    Nbps2 = Nbps/2;
    for kk = 1:K,
        tmp = data_in(:,kk).';
        
        % Slicing:
        tmps = 2*round((tmp-1-1i)/2)+1+1i;
        tmp_I = real(tmps);
        tmp_Q = imag(tmps);
        tmp_I(find(tmp_I>2^ceil(Nbps/2)-1)) = 2^ceil(Nbps/2)-1;
        tmp_I(find(tmp_I<-2^ceil(Nbps/2)+1)) = -2^ceil(Nbps/2)+1;
        tmp_Q(find(tmp_Q>2^floor(Nbps/2)-1)) = 2^floor(Nbps/2)-1;
        tmp_Q(find(tmp_Q<-2^floor(Nbps/2)+1)) = -2^floor(Nbps/2)+1;
        tmps = tmp_I+1i*tmp_Q;
        
        qi = real(tmp);
        qq = imag(tmp);
        qis = real(tmps);
        qqs = imag(tmps);
        y = [];
        
        %%% I component %%%
        for n = 1:Nbps2
            sb = (qi>0);
            ss = -2*sb+1;
            amp = (abs(qis)+1)/2;
            y = [y;amp.*(qi+ss.*(amp-1))];
            qi = ss.*(qi+(2^(Nbps2-n)).*ss);
            qis = ss.*(qis+(2^(Nbps2-n)).*ss);
        end;
        
        %%% Q component %%%
        for n = 1:Nbps2
            sb = (qq>0);
            ss = -2*sb+1;
            amp = (abs(qqs)+1)/2;
            y = [y;amp.*(qq+ss.*(amp-1))];
            qq = ss.*(qq+(2^(Nbps2-n)).*ss);
            qqs = ss.*(qqs+(2^(Nbps2-n)).*ss);
        end;
        
        data_out(:,kk) = y(:);
    end;
end;
data_out = sign(data_out);
end
