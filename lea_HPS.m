% function [f0_HPS, time_HPS, ampl_HPS] = lea_HPS(p,f,t)
% 
% This pitch tracking function is based on the Harmonic Product Spectrum
% (HPS) algorithm, as described un (Cuadra2001) using the spectrogram as 
% time-framing and DFT function.
% Need to be followed by function pitch_track_segments.m to convert detections
% into tracks.
%
    % INPUTS:
    %  - spectrogram p,f,t
    %
    % OUTPUTS:
    %  - f0: frequency estimate,
    %  - time: time vector of the frequency tracks,
    %  - ampl: amplitude of the tracks.

function [f0, time, ampl] = lea_HPS(p,f,t)

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
f0 = zeros(1,length(t));
ampl = zeros(1,length(t));

for i = 1:length(t)
    [M,I] = max(p_HPS(:,i));
    f0(i) = f_HPS(I);
    ampl(i) = sqrt(M);
end
time = t;
