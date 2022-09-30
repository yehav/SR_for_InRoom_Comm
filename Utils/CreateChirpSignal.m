function [chirp_signal, t] = CreateChirpSignal(fs, f_edges,...
                                               time_period_sec, total_signal_duration_sec)
% Input:
% * fs - sampling frequency [Hz]
% * f_edges - 1x2 vector - [f_min, f_max] [Hz]
% * time_period_sec - double
% * total_signal_duration_sec - double
% 
% Output:
% * chirp_signal - of length 'total_signal_duration_sec'
% * t - time vector [sec]
% 
% For Debug:
% fs = 16e3;% [Hz]
% f_edges = [200 7800];% [Hz]
% time_period_sec = 1;
% total_signal_duration_sec = 3;
% 
% [chirp_signal, t] = CreateChirpSignal(fs, f_edges, time_period_sec, total_signal_duration_sec)
% 
% Author: Yehav Alkaher.

t = 0:(1/fs):(total_signal_duration_sec - 1/fs);
phase = f_edges(1)*t + 0.5*t.^2*(diff(f_edges)/time_period_sec);
chirp_signal = db2mag(-40)*(1*sin(2*pi*phase)).*(t < time_period_sec) + 1*db2mag(-200)*randn(size(t));

end