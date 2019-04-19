% function [f0_HPS, time_HPS, ampl_HPS] = lea_HPS(p,f,t)
% 
% This pitch tracking function is based on the Harmonic Product Spectrum
% (HPS) algorithm, as described un (Cuadra2001) using the spectrogram as 
% time-framing and DFT function.
%
    % INPUTS:
    %  - spectrogram p,f,t
    %
    % OUTPUTS:
    %  - f0_HPS: Pitch track from HPS algorithm,
    %  - time_HPS: time of the HPS pitch track,
    %  - ampl_HPS: amplitude of the tracks.

function [f0_HPS, time_HPS, ampl_HPS] = lea_HPS(p,f,t)
% function [f0_HPS, time_HPS, ampl_HPS] = lea_HPS(x,fs,fft_size,overlap)
% [~,f,t,p] = spectrogram(x,hann(fft_size),round((overlap/100)*fft_size),fft_size,fs);

p_HPS = ones(length(f),length(t)); % initailization of p_HPS matrix to 1
nb_harmo = 2; % number of times the spctrum is downsampeled. In our case is to low to do it more thant twice. 
% For each time bin, the spectrum is downsampeled nb_harmo times. We keep
% the smallest frequency vector size. and complete the rest by nan values
% then the product of p_HPS(:,i) spectrum and the downsampled spectrum is
% realized.

for i = 1:length(t)
    h = p(:,i) ;% spectre au temps i
    for n = nb_harmo:-1:1
        h_down = downsample(h,n);
        diff_size = length(f)-length(h_down);
        h_down_sized = [h_down; NaN(diff_size,1)];
        p_HPS(:,i) = p_HPS(:,i).*h_down_sized;
    end
end
[a,~] = size(p_HPS);
f_HPS = linspace(0,50,a)*nb_harmo; %the associated frequency vector

% find Max frequency of the HPS, and the associated amplitude
f0_HPS = zeros(1,length(t));
ampl_HPS = zeros(1,length(t));

for i = 1:length(t)
    [M,I] = max(p_HPS(:,i));
    f0_HPS(i) = f_HPS(I);
    ampl_HPS(i) = sqrt(M);
end
time_HPS = t;