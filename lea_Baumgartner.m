% function [f0_baumgartner, time_baumgartner, ampl_baumgartner] = lea_Baumgartner(p,f,t,Ts)
% 
% This pitch tracking function is based on the algoritm described in
% Baumgartner2011. First a "thresholding is applied. Then, to find the
% pitch tracks, then a cost calculation is realized to associates STFT 
% pixels of a track together and have only one f0 value/time. 
%
    % INPUTS:
    %  - spectrogram p,f,t
    %  - Ts: fix threshold in dB.
    %
    % OUTPUTS:
    %  - f0_inst_freq: Pitch track from instantaneous frequency estimate,
    %  - time_inst_freq: time of the HPS pitch track,
    %  - ampl_inst_freq: amplitude of the tracks.

function [f0_baumgartner, time_baumgartner, ampl_baumgartner] = lea_Baumgartner(p,f,t,Ts)

threshold = 10.^(Ts/10); % real amplitude threshold

% Application of the threshold to the data
p_threshold = NaN(size(p));
for i = 1:length(f)
    ind = find(p(i,:)>=threshold);
    p_threshold(i,ind) = p(i,ind);
end

% ---------------------- COST CALCULATION -------------------- %
% The cost is calculated in two steps. First for the first "detected" value.
% Then for the following ones. % i is the time index, that starting at the first detected point (i0).
% <=> increase columns.

% Create vectors of indexes for  both time and frequency index for
% "detected" calls
[line,col] = find(isnan(p_threshold)==0); % line (= freq = j) and column (= time = i) index of non-NaN values
weight = 20 ; % User define weight

% Empty matrix declaration
Cost = zeros(1,length(t));
Detected_freq = NaN(1,length(t));
Detected_amp = NaN(1,length(t));
% For the first detected value 
% i = i0 + 1 
i = 2; j = 2; last_j = j;
P = weight * abs(log2(f(line(j-1))/f(line(j))));
Cost(i) = P - p_threshold(line(j),col(i));

% For the rest of the temporal observation
% i > i0 + 1 We go foward in time
for i = 3:length(t)-2
    % Are they any "detected frequency" ?
    frequency_index =  find(isnan(p_threshold(:,i))==0);
    % If they are, the following developpment aims to find the frequency
    % that minimize the cost, (especially in the case that there are
    % several simultaneous detected frequencies). 
    
    if isempty(frequency_index) == 0
        K = length(frequency_index); % Number of "Detected frequencies" 
        P = zeros(1,K); % P vector
        A = p_threshold(frequency_index,i)'; % Amplitudes
        for k = 1:K
            P(k) = weight * abs(log2(f(last_j)/f(frequency_index(k))));
        end
        Cost(i) = min((Cost(i-1) + P - A)); % minimum cost value
        if isnan(Cost(i))== 0
        ind_freq_Cost_min = find( (Cost(i-1) + P - A) == Cost(i),1); % Corresponding frequency
        Detected_freq(i) = f(frequency_index(ind_freq_Cost_min));
        Detected_amp(i) = A(ind_freq_Cost_min);
        last_j = ind_freq_Cost_min;
        end
    end
   
end 

for i = length(t)-2:-1:1
end

%  % ---------------------- -------------------- -------------------- %


f0_baumgartner = Detected_freq;
time_baumgartner = t;
ampl_baumgartner = Detected_amp;