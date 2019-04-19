% function [time_yin, f0_yin] = yin_estimator(p,f,t,x, fs,  varargin)

% Yin estimator function + segmentation of the results

function [time_yin, f0_yin, Start, End ] = lea_yin_estimator(x,fs,delta_t,delta_f,signal_mini_duration_s,varargin)
% from https://github.com/orchidas/Pitch-Tracking/blob/master/yin_estimator.m
% function that implements YIN algorithm for
% fundamental pitch tracking
% x - input audio signal
% fs - sampling rate
% time_yin,f0_yin - time_yin vector and associated fundamental frequencies estimated

% window size -  we assume the minimum f0_yin to be 1/0.25 = 4Hz 
win = round(0.25*fs);
N = length(x);
nframes = ceil(N/win);
% zero pad signal to have enough frames
[x] = vector_orientation(x,'column');
x = [x, zeros(1,win*nframes - N)];
x_frame = zeros(nframes, win);
start = 1;
%break into windows
for i = 1:nframes
    x_frame(i,:) = x(start:start + win - 1);
    start = start + win;
end

%step 1 - calculate difference function 
d = zeros(nframes,win);
x_temp = [x_frame, zeros(nframes,win)];
for tau = 0:win-1
    for j = 1:win  
         d(:,tau+1) = d(:,tau+1) + (x_temp(:,j) - x_temp(:,j+tau)).^2;         
    end
end

%step 2 - cumulative mean normalised difference function
d_norm = zeros(nframes,win);
d_norm(:,1) = 1;

for i = 1:nframes
    for tau = 1:win-1
        d_norm(i,tau+1) = d(i,tau+1)/((1/tau) * sum(d(i,1:tau+1)));
    end
end

% figure(1);
% subplot(211);
% plot(0:length(x)-1, reshape(d',1,length(x)));grid on;
% xlabel('Lags');
% ylabel('Difference function');
% subplot(212);
% plot(0:length(x)-1, reshape(d_norm',1,length(x)));grid on;
% xlabel('Lags');
% ylabel('Cumulative mean difference function');

%step 3 - absolute thresholding
lag = zeros(1,nframes);
th = 0.1;
for i = 1:nframes
    l = find(d_norm(i,:) < th,1);
    if(isempty(l) == 1)
        [v,l] = min(d_norm(i,:));
    end
    lag(i) = l; 
end

%step 4 - parabolic interpolation
period = zeros(1,nframes);
start = 1;

for i = 1:nframes
    if(lag(i) > 1 && lag(i) < win)
        alpha = d_norm(i,lag(i)-1);
        beta = d_norm(i,lag(i));
        gamma = d_norm(i,lag(i)+1);
        peak = 0.5*(alpha - gamma)/(alpha - 2*beta + gamma);
        %ordinate needs to be calculated from d and not d_norm - see paper
        %ordinate = d(i,lag(i)) - 0.25*(d(i,lag(i)-1) - d(i,lag(i)+1))*peak;
    else
        peak = 0;
    end
    %1 needs to be subtracted from 1 due to matlab's indexing nature
    period(i) = (lag(i)-1) + peak;
%     f0_yin(i,:) = fs/period(i)*ones(1,win);
%     time_yin(i,:) = ((i-1)*win:i*win-1)/fs;
end
f0_yin = fs./period;

time_yin_interm =linspace(0,length(x)/fs,length(f0_yin));
% time_yin_interm = time_yin_interm(f0_yin>10);
% f0_yin = f0_yin(f0_yin>10);
time_yin_interm(f0_yin<010) = NaN;
f0_yin(f0_yin<10) = NaN;


%for silent frames estimated frequency should be 0Hz
if ~isempty(varargin)
    [f0_yin] = silent_frame_classification(x_frame, f0_yin);
end

%step 5 - "Smooth" the instantaneous frequency using a median filter
f0_yin = medfilt1(f0_yin,11);

if isempty(f0_yin(isnan(f0_yin)==0))
    f0_yin = NaN;
    time_yin = NaN;
    Start = NaN;
    End = NaN;
else
%step 6 - separate the tracked tonals
Start = 1;
End = [];

% Detection of the pitch track contours (start and end point)
index = isnan(f0_yin)==0; % index = 0 <=> signal detecte
f0_short = f0_yin(index);
t_short = time_yin_interm(index);
for i = 2:length(f0_short)-1
    if (abs(f0_short(i)-f0_short(i-1))>delta_f)||(abs(t_short(i)-t_short(i-1))>delta_t),
        Start = [Start i];
        End = [End i-1];
    end
end
End = [End length(index)];

% Clear the detection of signal shorter than "signal_mini_duration"
for i = 1: length(Start)
    if (End(i)-Start(i) < round( signal_mini_duration_s/(time_yin_interm(2)-time_yin_interm(1)) ) )%|| (abs((f0_short(Start(i))-f0_short(End(i)))/(t_short(End(i))-t_short(Start(i))))<0.5),
        f0_short(Start(i):End(i))=NaN(1,End(i)-Start(i)+1);
        t_short(Start(i):End(i))=NaN(1,End(i)-Start(i)+1);
        End(i) = NaN;
        Start(i) = NaN;
    end
end
End = End(isnan(End)==0);
Start = Start(isnan(Start)==0);

Start_new = Start;
End_new = End;

% Attach segments that belong together
for i = 2:length(Start)
    if t_short(Start(i))-t_short(End(i-1)) < delta_t
        Start_new(i) = NaN;
        End_new(i-1) = NaN;
    end
end

End = End_new(isnan(End_new)==0);
Start = Start_new(isnan(Start_new)==0);
End(end) = length(t_short);

ind = find(End-Start<=0);
if isnan(ind)==0
    End = [End(1:ind-1) End(ind+1:end)];
    Start = [Start(1:ind-1) Start(ind+1:end)];
end

f0_yin = f0_short;
    % Good size
    time_yin = (0:N-1)/fs;
if length(Start) > 1

    f0_yin_final = NaN(size(time_yin));
    for i = 1:length(Start)
        Start_old = Start(i);
        End_old = End(i);
        Start(i) = find(time_yin >= t_short(Start_old),1);
        End(i) = find(time_yin >= t_short(End_old),1)-1;
        f0_yin_final(Start(i):End(i)) = interp1(t_short(Start_old:End_old),f0_yin(Start_old:End_old),time_yin(Start(i):End(i)));
    end
    f0_yin = f0_yin_final;

end
end



ind = find(End-Start<=0);
    if isnan(ind)==0
    End(ind) = NaN;
    Start(ind)= NaN;
    End = End(isnan(End)==0);
    Start = Start(isnan(Start)==0);
    end

