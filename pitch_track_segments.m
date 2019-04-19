function [f0,time,ampl,Start,End] = pitch_track_segments(f0,time,ampl,delta_t,delta_f,signal_mini_duration)
Start = 1;
End = [];

% delta_t = 3;%(s)
% delta_f = 0.8;%(Hz)
% signal_mini_duration = 15; %(samples)

% Detection of the pitch track contours (start and end point)
index = isnan(f0)==0; % index = 0 <=> signal detecte
f0_short = f0(index);
ampl_short = ampl(index);
t_short = time(index);
for i = 2:length(f0_short)-1
    if (abs(f0_short(i)-f0_short(i-1))>delta_f)||(abs(t_short(i)-t_short(i-1))>delta_t),
        Start = [Start i];
        End = [End i-1];
    end
end
End = [End length(index)];

% Clear the detection of signal shorter than "signal_mini_duration"
for i = 1: length(Start)
%    disp(i)
    if (End(i)-Start(i) < signal_mini_duration), %|| (abs((f0_short(Start(i))-f0_short(End(i)))/(t_short(End(i))-t_short(Start(i))))<0.5),
        f0_short(Start(i):End(i))=NaN(1,End(i)-Start(i)+1);
        ampl_short(Start(i):End(i))=NaN(1,End(i)-Start(i)+1);
        t_short(Start(i):End(i))=NaN(1,End(i)-Start(i)+1);
        End(i) = NaN;
        Start(i) = NaN;
    end
end
End = End(isnan(End)==0);
Start = Start(isnan(Start)==0);
f0 = f0_short;
ampl = ampl_short;
time = t_short;


% Attach segments that belong together
Start_new = Start;
End_new = End;

for i = 2:length(Start)
    if t_short(Start(i))-t_short(End(i-1)) < delta_t
        Start_new(i) = NaN;
        End_new(i-1) = NaN;
    end
end

End = End_new(isnan(End_new)==0);
Start = Start_new(isnan(Start_new)==0);
if isempty(End)==0
    End(end) = length(t_short);
end


