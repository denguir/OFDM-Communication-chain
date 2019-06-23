%% Lab 1:
% - Extract the impulse response at each position
% - Evaluate PDP and delay spread
% - Evaluate coherence bandwidth
% - Redo the exercise for BW=20MHz
clear variables;
clc;
close all;
% Load frequency response of 1000 local areas:
LOS = 1;
if LOS
    load('LOS.mat');
else
    load('NLOS.mat');  
end

% Baseband signal of BW=200MHz:
fs = 200e6; % Frequence sampling
ts = 1/fs; % Sampling period

%Baseband signal of BW=20MHz:
fs_low = 20e6;
ts_low = 1/fs_low;
ratio = fs/fs_low; % downsampling ratio
nbins = 501; % number of bins
fres = 4e5; % frequency resolution
fc = 2.35e9; % carrier frequency
bw = [fc-250*fres:fres:fc+250*fres]; % bandwidth
tau_bw = [0:nbins-1]*ts; % delay

% Evaluate impulse response of each local area:
tap_max = 6;
h_mat = zeros(10,10,10,ratio*tap_max);
h_mat_rect = zeros(10,10,10,tap_max);

for x=1:xx
    for y=1:yy
        for z=1:zz
            Hxyz = Data{1,x}{1,y}{1,z};
            hxyz = sqrt(nbins)*ifft(ifftshift(Hxyz));
            [~,start_index] = max(abs(hxyz));
            h_mat(x,y,z,:) = hxyz(start_index:start_index+ratio*tap_max-1);
            
            % downsampled version
            Hxyz_down = zeros(length(Hxyz), 1);
            Hxyz_down(225:275) = Hxyz(225:275);
            hxyz_down = sqrt(nbins)*ifft(ifftshift(Hxyz_down));
            [~,start_index] = max(abs(hxyz_down(1:end-20)));
            hxyz_down = hxyz_down(start_index:ratio:end);
            h_mat_rect(x,y,z,:) = hxyz_down(1:tap_max);
        end
    end
end

% Narrowband impulse response:
h_mat_rect_narrow = sum(h_mat_rect, 4);

% Evaluate PDP from impulse response:
% mean over all local area's impulse response
PDP = reshape(mean(mean(mean(abs(h_mat).^2, 1),2),3), [1, ratio*tap_max]);
PDP_rect = reshape(mean(mean(mean(abs(h_mat_rect).^2, 1),2),3), [1, tap_max]);

% Evaluate delay spread and coherence bandwidth:
tau = (0:tap_max*ratio-1)*ts;
tau_low = (0:tap_max-1)*ts_low;
[coherence_bw, sig_tau] = get_coh_bw(tau, PDP);
[coherence_bw_rect, sig_tau_rect] = get_coh_bw(tau_low, PDP_rect);
% fit PDP with exp
PDP_fit = 10*log10(PDP(1)*exp(-tau/sig_tau));
PDP_fit_rect = 10*log10(PDP_rect(1)*exp(-tau_low/sig_tau_rect));

figure;
plot(tau, 10*log10(PDP));
hold on; grid on;
plot(tau, PDP_fit);
xlabel('\tau [s]');
ylabel('PDP [dB]');
title('PDP @ 200MHz - PDP vs tau');
legend('Experimental PDP', 'Theoretical PDP');

figure;
plot(tau_low, 10*log10(PDP_rect));
hold on; grid on;
plot(tau_low, PDP_fit_rect);
xlabel('\tau [s]');
ylabel('PDP [dB]');
title('PDP @ 20MHz - PDP vs tau');
legend('Experimental PDP', 'Theoretical PDP');

figure;
plot(tau, abs(squeeze(h_mat(1,1,1,:))));
grid on;
xlabel('\tau [s]');
ylabel('|h|');
title('Impulse response at location (0,0,0)');


figure;
plot(tau_low, abs(squeeze(h_mat_rect(1,1,1,:))));
grid on;
xlabel('\tau [s]');
ylabel('|h|');
title('Impulse response at location (0,0,0)');

% Show original transfer function
figure;
plot(bw, 20*log10(abs(Data{1,1}{1,1}{1,1})));
grid on;
xlabel('f [Hz]');
ylabel('|H| [dB]');
title('Transfer function at location (0,0,0)');

% Show original impulse response 
figure;
plot(tau_bw, abs(sqrt(nbins)*ifft(ifftshift(Data{1,1}{1,1}{1,1}))));
grid on;
xlabel('\tau [s]');
ylabel('|h|');
title('Impulse response at location (0,0,0)');

%% Lab 2:
% - Construct narrowband channel + describe its distribution
% - In wideband (20MHz), describe distribution of each tap 
% - Build statistical channel model for narrowband and wideband

% Narrowband channel:
% one tap per impulse response:
narrow_taps = abs(h_mat_rect_narrow);
narrow_taps = reshape(narrow_taps, [], 1);
% Narrowband statistical distribution (LOS & NLOS):
if LOS
    % - LOS (Rician):
    rice_narrow = fitdist(narrow_taps, 'Rician');
    K_rice_narrow = rice_narrow.s^2/(2*rice_narrow.sigma^2);
else
    % - NLOS (Rayleigh):
    ray_narrow = fitdist(narrow_taps, 'Rician');
    K_ray_narrow = ray_narrow.s^2/(2*ray_narrow.sigma^2);
%     ray_narrow = fitdist(narrow_taps, 'Rayleigh');
%     K_ray_narrow = ray_narrow.B;
end

% Wideband channel:
wide_taps = abs(h_mat_rect);
wide_taps = reshape(wide_taps, [], tap_max);
% Wideband statistical distribution (LOS & NLOS):
if LOS
    % - LOS (Rician)
    rice_wide = cell(1, tap_max);
    K_rice_wide = zeros(1, tap_max);
    for i=1:tap_max
        if i < 3
            rice_wide{i} = fitdist(wide_taps(:,i), 'Rician');
            K_rice_wide(i) = rice_wide{i}.s^2/(2*rice_wide{i}.sigma^2);
        else
            rice_wide{i} = fitdist(wide_taps(:,i), 'Rayleigh');
            K_rice_wide(i) = rice_wide{i}.B;
        end
    end
else
    % - NLOS (Rayleigh):
    ray_wide = cell(1, tap_max);
    K_ray_wide = zeros(1, tap_max);
    for i=1:tap_max
          ray_wide{i} = fitdist(wide_taps(:,i), 'Rician');
          K_ray_wide(i) = ray_wide{i}.s^2/(2*ray_wide{i}.sigma^2);
%         ray_wide{i} = fitdist(wide_taps(:,i), 'Rayleigh');
%         K_ray_wide(i)= ray_wide{i}.B;
    end
end

%% Lab 3:
% - Implement parr transceiver
% - Assess BER of the channel modeled only by AWGN
% - Assess BER of the channel taking into account the channel impulse
% response
Nbits = 2^15; % number of bits to transmit
Nchannels = 64; % number of OFDM sub-channels
size_cp = 16; % cycle prefix size, avoids IBI
modulation = 'qam';
Nbps = 4; % number of bits per symbol
%noise_seed = rng(0); % seed used for the noise to compare BER curve 

% Transmission of bits
bits = randi([0, 1], Nbits, 1); % random bits (0, 1)
% insert pilot bits for CFO detection
pilot_idx = [-21, -7, 7, 21] + 32; % [-21, -7, 7, 21] index used in Wifi
bits = insert_pilot(bits, pilot_idx, Nbps, Nchannels);
bitsPhase = bits;
bitsPhase(bitsPhase==0) = -1; % random bits (-1, 1)

% Bits to Symbols (in freq domain):
symb = mapping(bits, Nbps, modulation);
symb_pilot = symb(pilot_idx);
parr_symb = reshape(symb, Nchannels, []); % parallelize parr symbols
time_symb = ifft(parr_symb, [], 1); % Get symbols in time domain

% Add cyclic prefix to each block of symbol:
prefixed_symb = repmat(time_symb, 2, 1);
prefixed_symb = prefixed_symb(Nchannels - size_cp + 1:end, :);
serial_symb = reshape(prefixed_symb, [], 1);

NR = 10; % number of realisations
% We tweak the value of the SNR (EbN0) to build BER curves
EbN0 = -5:2:30; % in dB
BER = zeros(NR, length(EbN0));
%%
for rr=1:NR % run for NR channel realisations
    % Define the groundtruth wideband channel impulse response:
    if LOS
        % LOS:
        h_wide = get_ch_rsp(rice_wide);
        %h_wide = [1 0 0 0 0 0]; % AWGN
        H_wide = fft(h_wide, Nchannels).'; % pad with zeros h_wide_los
    else
        %NLOS:
        h_wide = get_ch_rsp(ray_wide);
        H_wide = fft(h_wide, Nchannels).'; 
    end

    % Convolve transmitted signal with the channel
    % apply channel impulse response:
    % and add AWGN to the transmitted signal:
    for k=1:length(EbN0)
        r = conv(serial_symb, h_wide);
        r = r(1:end-length(h_wide)+1); % remove tail of convolution
        % add noise:
        signal_power = sum(abs(r).^2)*ts_low;
        Eb = 0.5 * signal_power/Nbits;
        N0 = Eb/(10.^(EbN0(k)/10));
        noise_power = 2*N0*fs_low;
        %rng(noise_seed);
        noise = sqrt(noise_power/2)*((randn(length(r),1))+1i*(randn(length(r),1)));    
        r = r + noise;
        
        % parallelize:
        parr_rcv_symb = reshape(r, Nchannels + size_cp, []);
        rcv_symb = parr_rcv_symb(size_cp+1:end, :); % remove cycle prefix
        % switch in frequency domain for comparison with sent signal:
        freq_rcv_symb = fft(rcv_symb, [], 1);
        % equalization of each sub-carrier independently:
        equalized = freq_rcv_symb./H_wide;
        serial_rcv_symb = reshape(equalized, [], 1);
        bitsEstimated = demapping(serial_rcv_symb, Nbps);
        BER(rr, k) = length(find(bitsPhase~=bitsEstimated))/length(bits);
    end
end
figure;
BER_final=mean(BER, 1);
semilogy(EbN0, BER_final);
title('BER curve, channel known');
grid on; 

%% Lab 4:
% Estimate the channel by using known preamble symbols
% Evaluate channel bitsEstimated wrt real channel impulse response

% Define the preamble symbols, used for channel bitsEstimated
% Its 2 times the same preamble (useful for ToA bitsEstimated)
power_symb = sum(abs(symb).^2)/length(symb);
preamble_symb = randi([0, 1], Nchannels, 1);
preamble_symb(preamble_symb==0) = -1;
preamble_symb = preamble_symb * sqrt(power_symb);
preamble_symb_time = ifft(preamble_symb);

preamble_tot = repmat(preamble_symb_time, [2, 1]);
preamble_tot = [preamble_tot(end - 2 * size_cp + 1:end); preamble_tot];
% construction of frame (preamble + data with cycle prefix):
frame_symb = [preamble_tot; serial_symb];

% NMSE is the vector which evaluates the error between the estimated
% channel and the true channel for each EbN0
NMSE_H = zeros(size(BER));
%%
for rr=1:NR
    % Define the groundtruth wideband channel impulse response:
    if LOS
        % LOS:
        h_wide = get_ch_rsp(rice_wide);
        %h_wide = [1 0 0 0 0 0];
        H_wide = fft(h_wide, Nchannels).'; % pad with zeros h_wide_los
    else
        %NLOS:
        h_wide = get_ch_rsp(ray_wide);
        %h_wide = [1 0 0 0 0 0];
        H_wide = fft(h_wide, Nchannels).'; 
    end
    for k=1:length(EbN0)
        % apply channel impulse response:
        r = conv(frame_symb, h_wide);
        r = r(1:end-length(h_wide)+1); % remove tail of convolution
        % add noise:
        signal_power = sum(abs(r).^2)*ts_low;
        Eb = 0.5 * signal_power/Nbits;
        N0 = Eb/(10.^(EbN0(k)/10));
        noise_power = 2*N0*fs_low;
        % rng(noise_seed);
        noise = sqrt(noise_power/2)*((randn(length(r),1))+1i*(randn(length(r),1)));    
        r = r + noise;

        % parallelize:
        parr_rcv_symb = reshape(r, Nchannels + size_cp, []);
        % the symbols corresponding to the preamble are contained in the two
        % first columns of parr_rcv_symb, because: lenght(2 * CP + 2 *
        % preamble) = 160 and one column has lenght size_cp + N_channels = 80:
        parr_preamble_rcv = parr_rcv_symb(:, 1:2);
        serial_preamble = reshape(parr_preamble_rcv, [], 1); % serialize
        serial_preamble = serial_preamble(2*size_cp+1:end); % remove cycle prefix
        parr_preamble_rcv = reshape(serial_preamble, Nchannels, []);

        parr_preamble_freq = fft(parr_preamble_rcv, [], 1);
        H_wide_est = parr_preamble_freq ./ preamble_symb;
        H_wide_est = mean(H_wide_est, 2);
        h_wide_est = ifft(H_wide_est);
        h_wide_est = h_wide_est(1:6);

        rcv_symb = parr_rcv_symb(size_cp+1:end, 3:end); % remove cycle prefixed
        % switch in frequency domain for comparison with sent signal:
        freq_rcv_symb = fft(rcv_symb, [], 1);
        % equalization of each sub-carrier independently:
        equalized = freq_rcv_symb./H_wide_est;
        serial_rcv_symb = reshape(equalized, [], 1);
        bitsEstimated = demapping(serial_rcv_symb, Nbps);
        NMSE_H(rr, k) = sum(abs(H_wide_est-H_wide).^2)/sum(abs(H_wide).^2);
        BER(rr, k) = length(find(bitsPhase~=bitsEstimated))/length(bits);
    end
end
figure;
BER1 = mean(BER, 1);
semilogy(EbN0,BER1); 
title('BER curve, channel estimated');
grid on;

figure;
NMSE_H1 = mean(NMSE_H, 1);
semilogy(EbN0,NMSE_H1);
title('NMSE between H_{est} and H_{true}');
grid on;

%% Lab 5 & Lab 6:
% Add uncertainty on time of arrival (ToA)
% Explain the phase rotation in H implied by ToA
% Implement time acquisition to estimate ToA
% Evaluate the accuracy of ToA bitsEstimated (BER performance)
% Add CFO to the communication
% Compensation and tracking of CFO
% Evaluate BER and Constellation with CFO

ToA = 20; % 30 symbols = (30 * ts) seconds
CFO = 0*fs_low*1e-6; % in ppm
% delayed_frame_symb = circshift(frame_symb, ToA);
errorToA = zeros(size(BER));
RMSEToA = zeros(1, size(BER, 2));

for rr=1:NR
    % Define the groundtruth wideband channel impulse response:
    if LOS
        % LOS:
        h_wide = get_ch_rsp(rice_wide);
        %h_wide = [1 0 0 0 0 0];
        H_wide = fft(h_wide, Nchannels).'; % pad with zeros h_wide_los
    else
        %NLOS:
        h_wide = get_ch_rsp(ray_wide);
        %h_wide = [1 0 0 0 0 0];
        H_wide = fft(h_wide, Nchannels).';
    end

    for k=1:length(EbN0)
        % apply channel impulse response:
        r = conv(frame_symb, h_wide);
        r = r(1:end-length(h_wide)+1); % remove tail of convolution
        % add noise:
        signal_power = trapz(abs(r).^2)*ts_low;
        Eb = 0.5 * signal_power/Nbits;
        N0 = Eb/(10.^(EbN0(k)/10));
        noise_power = 2*N0*fs_low;
        % rng(noise_seed);
        noise = sqrt(noise_power/2)*((randn(length(r),1))+1i*(randn(length(r),1)));  
        r = r + noise;
       
        r = circshift(r, ToA);
        % cfo:
        phase_cfo = 2*pi*CFO*ts_low*(0:length(r)-1)';
        r = r .* exp(-1j*phase_cfo); % add cfo to reveived signal

        % estimate ToA by analyzing correlation in the frame
        ToA_est = find_start_frame(r, Nchannels, 2*size_cp);
        errorToA(rr, k) = (ToA - ToA_est)^2;
        % correct ToA:
        r = circshift(r, -ToA_est+size_cp/2);
        % correct CFO:
%         CFO_est = find_cfo(r, Nchannels, 2*size_cp, ts_low);
%         phase_cfo_est = 2*pi*CFO_est*ts_low*(0:length(r)-1)';
%         r = r .* exp(1j*phase_cfo_est); % add cfo to reveived signal
        % parallelize:
        parr_rcv_symb = reshape(r, Nchannels + size_cp, []);
        % the symbols corresponding to the preamble are contained in the two
        % first columns of parr_rcv_symb, because: lenght(CP + 2 *
        % preamble) = 160 and one column has lenght size_cp + N_channels = 80:
        parr_preamble_rcv = parr_rcv_symb(:, 1:2);
        serial_preamble = reshape(parr_preamble_rcv, [], 1); % serialize
        serial_preamble = serial_preamble(2*size_cp+1:end); % remove cycle prefix
        parr_preamble_rcv = reshape(serial_preamble, Nchannels, []);

        parr_preamble_freq = fft(parr_preamble_rcv, [], 1);
        H_wide_est = parr_preamble_freq ./ preamble_symb;
        H_wide_est = mean(H_wide_est, 2);
        h_wide_est = ifft(H_wide_est);
        h_wide_est = h_wide_est(1:tap_max);
        
        rcv_symb = parr_rcv_symb(size_cp+1:end, 3:end); % remove cycle prefix
        
        % switch in frequency domain for comparison with sent signal:
        freq_rcv_symb = fft(rcv_symb, [], 1);
        
        % equalization of each sub-carrier independently:
        equalized = freq_rcv_symb./H_wide_est;    
        
        % frequency tracking (CFO)
%         phaseShift = mean(angle(equalized(pilot_idx, :)) - angle(symb_pilot), 1);
%         equalized = equalized .* exp(-1j*phaseShift);
        
        serial_rcv_symb = reshape(equalized, [], 1);
        bitsEstimated = demapping(serial_rcv_symb, Nbps);
        NMSE_H(rr, k) = immse(abs(H_wide), abs(H_wide_est))/sum(abs(H_wide).^2);
        BER(rr, k) = length(find(bitsPhase~=bitsEstimated))/length(bits);
    end
end

figure;
BER2=mean(BER, 1);
semilogy(EbN0,BER2);
title('BER curve, channel and ToA estimated');
grid on;

figure;
semilogy(EbN0,mean(NMSE_H, 1));
title('NMSE between H_{est} and H_{true} when ToA estimated');
grid on;

figure;
plot(real(reshape(freq_rcv_symb, [], 1)), imag(reshape(freq_rcv_symb, [], 1)), 'r*'); hold on;
plot(real(symb), imag(symb), 'b*');

figure;
NMSE_ToA = sqrt(mean(errorToA, 1));
plot(EbN0, NMSE_ToA);
title('RMSE ToA');
grid on;


%% Lab 7:
% Estimate main angles of arrival (AoA) of waves
% Use duality AoA - Spatial correlation to compute spatial correlation
% along directions X, Y and Z for narrowband and wideband models
% Compare the correlation with Clarke's Model

addpath(genpath('Lab7'));
c=3e+8; % speed of light
fc=2.35e+9; % carrier frequency
lambda=c/fc; % wave length
beta = 2*pi/lambda; % wave number
dist = 2*1e-2; % distance between each antenna

% Range of angles (spherical coordinates):
phi = (-1:0.025:1) * pi;
theta = (0:0.0125:1) * pi;

% choose between wideband and narrowband models
wide = 1;
if wide
    h = h_mat_rect;
else
    h = h_mat_rect_narrow;
end
% matrix of AoA:
AoA = angle_of_arrival(h, phi, theta, beta, dist);
maxAoA = 10*log10(2*pi*max(max(abs(AoA(:,:,1)).^2)));
minAoA = 10*log10(2*pi*min(min(abs(AoA(:,:,1)).^2)));
% Plot angular spectrum vs (phi, theta) 
% a0= a(21:-1:1,:,1); % phi from 0 to -pi(=pi)
% a1 = a(end:-1:22,:,1); % phi from pi to 07
Sn=10*log10(2*pi*abs(AoA(:,:,1)').^2);
plotImage(phi, theta, Sn , [maxAoA, minAoA]);
index_simo={};

for i=1:6
    if i <=2
        index_simo{1,i}=findLocalMaxima(10*log10(2*pi*abs(AoA(:,:,i)').^2),-70);  
    else
        index_simo{1,i}=findLocalMaxima(10*log10(2*pi*abs(AoA(:,:,i)').^2),-90);  
    end
end
%channel estimation
plot(phi(index_simo{1,1}(2,:)), theta(index_simo{1,1}(1,:)), 'o', 'LineWidth', 2);

% Channel estimation from main AoA, taking into account spatial phase shift
h_simo = get_simo_ch(AoA, index_simo, phi, theta, h, beta, dist);

% Evaluate spatial correlation:
dX = (0:0.1:xx) * dist;
dY = (0:0.1:yy) * dist;
dZ = (0:0.1:zz) * dist;
tap = 1;

main_phi_idx = index_simo{1,tap}(2,:);
main_theta_idx = index_simo{1,tap}(1,:);
beta_vec = beta*[sin(theta(main_theta_idx)).*cos(phi(main_phi_idx));
                      sin(theta(main_theta_idx)).*sin(phi(main_phi_idx));
                      cos(theta(main_theta_idx))];

    
main_AoA = diag(AoA(main_phi_idx, main_theta_idx , tap)); % main values of AoA
clarke = besselj(0, dZ*beta);
[Rx, Ry, Rz] = spatial_corr(main_AoA, beta_vec, dX, dY, dZ);

% normalize correlation and plot:
Rx_abs = abs(Rx)/abs(Rx(1));
Ry_abs = abs(Ry)/abs(Ry(1));
Rz_abs = abs(Rz)/abs(Rz(1));

figure;
plot(dZ, Rz_abs); hold on;
plot(dX, Rx_abs); hold on;
plot(dY, Ry_abs); grid on;
% plot(dZ, clarke); hold on;
%% lab 9:
% Implement SIMO communication with Maximum Ratio Combining (MRC)
% Evaluate the impact of the number of antennas and the distance between them

R_considered=[]; % matrix of positions
% place your antenna's wisely according to spatial correlation figure
R_considered(1,:)=[0 0 0]*dist;
R_considered(2,:)=[1 1 1]*dist;
% R_considered(2,:)=[8 4 9]*dist;

% R_considered(2,:)=[6 5 4]*dist;
% R_considered(3,:)=[4 8 2]*dist;


nb_antennas=size(R_considered,1);
if wide
    h_considered = zeros(nb_antennas, tap_max);
else
    h_considered = zeros(nb_antennas, 1);
end

H_est = zeros(nb_antennas, Nchannels);
Rsymb = zeros(nb_antennas, Nchannels, Nbits/(Nbps*Nchannels));
BER = zeros(NR,length(EbN0));

ToA = 0; % 30 symbols = (30 * ts) seconds
CFO = 0*fs_low*1e-6; % in ppm

for it=1:NR
    for nb=1:nb_antennas
        for j=1:size(h_considered,2)
            R = R_considered(nb,:);
            for i=1:size(index_simo{1,j},2)
                idx1=index_simo{1,j}(1,i); % list of phi idx
                idx2=index_simo{1,j}(2,i); % list of theta idx
                beta = (2*pi/lambda)*[sin(theta(idx2))*cos(phi(idx1)),sin(theta(idx2))*sin(phi(idx1)),cos(theta(idx2))];
                angle_rand=random('unif',0,2*pi);
                expo=exp(1i*(angle_rand-dot(beta,R)));
                h_considered(nb,j)=h_considered(nb,j)+AoA(idx1,idx2,j)*expo;
            end  
        end
    end
    %noise_base=(randn(length(frame_symb)+length(h_considered(1,:))-1,1))+1i*(randn(length(frame_symb)+length(h_considered(1,:))-1,1));
    for k=1:length(EbN0)
        for i=1:nb_antennas
            % Define the groundtruth wideband channel impulse response:
                % LOS:
            h_wide = h_considered(i,:);
            H_wide = fft(h_wide, Nchannels).'; % pad with zeros h_wide_los

            % apply channel impulse response:
            r = conv(frame_symb, h_wide);
            % add noise:
            signal_power = sum(abs(r).^2)*ts_low;
            Eb = 0.5 * signal_power/Nbits;
            N0 = Eb/(10.^(EbN0(k)/10));
            noise_power = 2*N0*fs_low;
            
            noise = sqrt(noise_power/2)*((randn(length(r),1))+1i*(randn(length(r),1)));    
            r = r + noise;
            r = r(1:end-length(h_wide)+1); % remove tail of convolution
            r = circshift(r, ToA);
            % cfo:
%             phase_cfo = 2*pi*CFO*ts_low*(0:length(r)-1)';
%             r = r .* exp(1j*phase_cfo); % add cfo to reveived signal

            % estimate ToA by analyzing correlation in the frame
%             ToA_est = find_start_frame(r, Nchannels, 2*size_cp);
%             % correct ToA:
%             r = circshift(r, -ToA_est+size_cp/2);
            % correct CFO:
%             CFO_est = find_cfo(r, Nchannels, 2*size_cp, ts_low);
%             phase_cfo_est = 2*pi*CFO_est*ts_low*(0:length(r)-1)';
%             r = r .* exp(-1j*phase_cfo_est); % add cfo to reveived signal
            % parallelize:
            parr_rcv_symb = reshape(r, Nchannels + size_cp, []);
            % the symbols corresponding to the preamble are contained in the two
            % first columns of parr_rcv_symb, because: lenght(CP + 2 *
            % preamble) = 160 and one column has lenght size_cp + N_channels = 80:
            parr_preamble_rcv = parr_rcv_symb(:, 1:2);
            serial_preamble = reshape(parr_preamble_rcv, [], 1); % serialize
            serial_preamble = serial_preamble(2*size_cp+1:end); % remove cycle prefix
            parr_preamble_rcv = reshape(serial_preamble, Nchannels, []);

            parr_preamble_freq = fft(parr_preamble_rcv, [], 1);
            H_wide_est = parr_preamble_freq ./ preamble_symb;
            H_wide_est = mean(H_wide_est, 2);
            h_wide_est = ifft(H_wide_est);
            h_wide_est = h_wide_est(1:6);

            rcv_symb = parr_rcv_symb(size_cp+1:end, 3:end); % remove cycle prefix

            % switch in frequency domain for comparison with sent signal:
            freq_rcv_symb= fft(rcv_symb, [], 1);
            
            H_est(i,:) = H_wide_est;
            Rsymb(i,:,:) = freq_rcv_symb;
            % equalization of each sub-carrier independently:
            %equalized = freq_rcv_symb./H_wide_est;
            %equalized = equalized + freq_rcv_symb.*conj(H_wide_est);    
            %normH = normH + abs(H_wide_est).^2;
        end
        % equalization of each sub-carrier independently using multiple antenna:
        norm_H = sum(abs(H_est).^2, 1)';
        equalized = zeros(Nchannels, Nbits/(Nbps*Nchannels));
        for n=1:nb_antennas
            equalized = equalized + squeeze(Rsymb(n,:,:)) .* H_est(n,:)';  
        end
        equalized = equalized ./ norm_H;
        % frequency tracking (CFO)
%         phaseShift = mean(angle(equalized(pilot_idx, :)) - angle(symb_pilot), 1);
%         equalized = equalized .* exp(-1j*phaseShift);
        
        serial_rcv_symb = reshape(equalized, [], 1);
        bitsEstimated = demapping(serial_rcv_symb, Nbps);
        BER(it,k) = length(find(bitsPhase~=bitsEstimated))/length(bits);
    end
end

figure;
BER_fin=mean(BER, 1);
semilogy(EbN0,BER_fin);
title('BER curve, SIMO');
grid on;
