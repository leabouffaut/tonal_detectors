% function [f0_inst_freq, time_inst_freq, ampl_inst_freq] = lea_inst_freq(p,f,t,x,fs)
% 
% This pitch tracking function is based on instantaneous frequency
% estimate in time. A smoothing median filter is applied as well as
% mean variance measurement for rejecting part of the data.
%
    % INPUTS:
    %  - spectrogram p,f,t
    %  - x: observation signal,
    %  - fs:  sampling frequency (Hz),
    %
    % OUTPUTS:
    %  - f0_inst_freq: Pitch track from instantaneous frequency estimate,
    %  - time_inst_freq: time of the HPS pitch track,
    %  - ampl_inst_freq: amplitude of the tracks.


function [f0_inst_freq, time_inst_freq, ampl_inst_freq] = lea_inst_freq(p,f,t,x,fs)
% function [f0_inst_freq, time_inst_freq, ampl_inst_freq] = lea_inst_freq(x,fs,fft_size,overlap)
%[~,f,t,p] = spectrogram(x,hann(fft_size),round((overlap/100)*fft_size),fft_size,fs);

% Instantaneous frequency measurement
z = hilbert(x);
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
[x] = vector_orientation(x,'line');

% "Smooth" the instantaneous frequency using a median filter
instfreq_filt = medfilt1(instfreq,21);

% Amplitude of these frequencies on the spectrogram
tx = (0:length(x)-1)/fs;
p_interp = interp2(t,f,p,tx,f);
ampl_inst_freq = NaN(size(instfreq_filt));
for i = 1:length(instfreq_filt)
    int_f = find(f>= instfreq_filt(i),1);
    if isempty(int_f)==0,
    ampl_inst_freq(i) = p_interp(int_f,i);
    end
end

% win_size = 501; % 5s
% var_f0 = zeros(size(instfreq));
% instfreq_ZP = [zeros(1,floor(win_size/2))  instfreq' zeros(1,floor(win_size/2))];
% for i = 1:length(var_f0),
%     temp = instfreq_ZP (i:i+win_size-1);
%    var_f0(i) = var(temp(isnan(temp)==0));
% end
% 
% f0_inst_freq = instfreq_filt(var_f0<50);
% time_inst_freq = tx(var_f0<50);
% ampl_inst_freq = ampl_inst_freq(var_f0<50);

f0_inst_freq = [instfreq_filt 0];
time_inst_freq = tx;
ampl_inst_freq = [ampl_inst_freq 0];

