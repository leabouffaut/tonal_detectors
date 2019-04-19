% x = generate_observation_SNR_controlled(noise,signal,fs,fsig,time_of_arrival,SNR)
%
% This fonction can be used to generate an observation with a controlled
% SNR, following this definition SNR = 10*log10(signal_power/noise_power)
% with x_power = x_energy/length(x). This fonction is build for tonals.
% Therefore, the energy term is relative to the energy in the signal's
% frequency band, define by fsig = [fmin fmax]. The noise energy is 
% measured at each time_of_arrival value for the duration of the signal,
% and the signal's amplitude is ajusted by the variable coef to satisfy the
% imposed SNR. This function is based on the method proposed (Mellinger & 
% Clark 2006, see Mobysound website: http://www.mobysound.org/software.html)
%
% INPUTS:
%    - noise, noise vector where the signal is going to be injected
%    with a controlled SNR, sampled at fs;
%    - signal, the signal that is going to be injected at each time_of_arrival
%    sampled at fs;
%    - fs, the sampling frequency (Hz);
%    - fsig, fsig = [fmin fmax], the signal's frequency band;
%    - time_of_arrival, time vector containig the starting time where
%    signals will be injected;
%    - SNR, imposed Signal to Noise Ratio  (dB);
%
% OUTPUT: 
%    - x, temporal observation vector of a mixture of signal and noise at
%    the required SNR
%    - real_SNR (optional), double-checked SNR value according 


function [x,real_SNR] = generate_observation_SNR_controlled(noise,signal,fs,fsig,time_of_arrival,SNR)
% pre-processing of the data
noise = noise - mean(noise); % so the noise is centered
signal = signal - mean(signal); signal = signal/max(abs(signal)); % the signal is normalized

%% Signal and noise mixture
% Initialize the received signal with the noise
received_signal = noise;

% We want to measure the exact noise power at the TOA where we're going to
% inject the signal, to adjust its power and have the required SNR.
N = length(signal);
N_2 = 2 .^ nextpow2(N);

% Frequency indexes corresponding to fmin and fmax of the simulated Z-call
% +/- 1 Hz (upper and lower)
i0 = round((fsig(1)-1)/(fs/2)*N_2/2) + 1;	% index for lower freq bound
i1 = round((fsig(2)+1)/(fs/2)*N_2/2) + 1;	% ...and upper

% Signal Energy
Ref_Signal_energy = energy_measurement(signal,[i0 i1],[N N_2]);

% Initialize 
real_SNR = [];
% For each emission
for emmission_number = 1: length(time_of_arrival)
    
    % Sample where we'll start injecting the call
    arrival_sample = round(time_of_arrival(emmission_number)*fs);
    % Noise cut from the arrival sample and for the duration of the call
    noise_short = received_signal(arrival_sample:arrival_sample+N-1);
    % Calculation of noise_short power in the bandwidth of the Z-call
    Noise_energy = energy_measurement(noise_short,[i0 i1],[N N_2]);
    % Calculation of the coefficient that will be used to set signal at
    % the right amplitude
    Signal_wanted_energy = Noise_energy * 10^(SNR/10) - Noise_energy;
    coef = sqrt(Signal_wanted_energy/Ref_Signal_energy);

    % Received signal
    received_signal(arrival_sample:arrival_sample+N-1) = received_signal(arrival_sample:arrival_sample+N-1) + coef*signal;
    
    if nargout == 2
    % SNR (to verify)
    Observation_energy = energy_measurement(received_signal(arrival_sample:arrival_sample+N-1),[i0 i1],[N N_2]);
    real_SNR = [real_SNR 10*log10(Observation_energy/Noise_energy)];
    end
    
end
% Returned built simulated signal
x = received_signal;
x = x-mean(x); x = x/max(x);

function energy = energy_measurement(x,I,NN)
pad = [x zeros(1, NN(2)-NN(1))]; % zero-pad to power-of-2 length
FT = abs(fft(pad));
energy = (1/NN(2)) * sum(FT(I(1):I(2)).^2)/NN(1);