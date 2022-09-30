%% Synopsis:
% * This script demonstrates
%   the temporal howling detection algorithm for Speech-Reinforcement Systems.
% * The closed-loop transfer-function (TF) of the SR system is simulated
%   via a two-pole TF.
%       
% * The script includes:
%   - Generate Sample Signals
%   - TF Design by Dominant Pole Configuration
%   - Apply TF on Test Signal
%   - Retrospective Howling Detection
%   - Analysis - Finding First Howling Detection
% 
% Author: Yehav Alkaher.

clc; close all; clear;
%% Initialization:
format shortg
addpaths_script

% Initialization Script:
Initialize_Speech_Reinforcement_System;

rng(noise_rng_seed);

% targetFolder = 'Results\HowlingTests\';
% if ~isdir(targetFolder)
%     mkdir(targetFolder)
% end

[MSD, GC] = GenerateMsdGainControlStructs(fs, 1);
%% Generate Sample Signals:
% Chirp:
total_signal_duration_sec = 3;
f_edges = [200 7800];
time_period_sec = 1;
[chirp_signal, t] = CreateChirpSignal(fs, f_edges, time_period_sec, total_signal_duration_sec);

% Low-level Gaussian Noise
noise_duration_sec = (MSD.Const.msd_frame_length + MSD.Const.msd_frame_shift)/fs;
low_level_noise_signal = [...
    db2mag(-40)*randn(1,noise_duration_sec*fs),...
    db2mag(-200)*randn(1,(total_signal_duration_sec - noise_duration_sec)*fs)];

% Speech:
Speech_dir_path = '';
Speech_filename = '564352__anzbot__ladies-and-gentlemen.mp3';
[y,Fs] = audioread([Speech_dir_path Speech_filename]);
y = y(:,1).';
if (fs ~= Fs)
    [U,D] = rat(fs/Fs);
    test_speech_signal = resample(y,U,D);
else
    test_speech_signal = y;
end
noiseStd = db2mag(-60)*std(test_speech_signal);
test_speech_signal = test_speech_signal + noiseStd*randn(size(test_speech_signal));
%% TF Design by Dominant Pole Configuration:
% Input:
% - sampling frequency
% - desired pole frequency
% - desired frequency-magnitude change-rate (increasing/decreasing)
% Output:
% - pole 
howl_freq_vec = 2000:250:6800;% Hz
howl_freq_vec = unique(round(howl_freq_vec/(fs/MSD.Const.msd_frame_length))*(fs/MSD.Const.msd_frame_length));

decreasing_howling_rates_dB_sec_vec = [-1000:100:-600, -500:50:-50];
decreasing_howling_pole_magnitudes_vec = db2mag(decreasing_howling_rates_dB_sec_vec/fs);

increasing_howling_rates_dB_sec_vec = [0:10:150, 200:100:1500];
increasing_howling_pole_magnitudes_vec = db2mag(increasing_howling_rates_dB_sec_vec/fs);

% zeta = decreasing_howling_pole_magnitudes_vec(10);
% zeta = increasing_howling_pole_magnitudes_vec(10);
% f_pole = howl_freq_vec(1);
zeta = db2mag(-300/fs);
f_pole = 2000;% Hz

% Chirp-wise:
[~, chirp_at_f_pole_index] = min(abs(f_edges(1) + t(t < time_period_sec)*(diff(f_edges)/time_period_sec) - f_pole));

z_pole = exp(1i*2*pi*f_pole/fs);
P_2 = tf(1, real([1, -zeta*(z_pole + 1/z_pole), zeta^2*1]), 1/fs, 'variable','z^-1');
fvtool(P_2.Numerator{1}, P_2.Denominator{1}, 'polezero')

assert(all(...
         abs(abs(roots(P_2.Denominator{1})) - zeta)...
         < 1e-5))

%% Apply TF on Test Signal:
signal = [...
    1*db2mag(-200)*randn(1,1.5*fs),...
    low_level_noise_signal,...
    chirp_signal...
    test_speech_signal...
    ];
noisy_signal = signal + 0*db2mag(-150)*randn(size(signal));

sample_index_to_filter_index_mat = [[1/fs, 1.5, 3, 4]*fs;[2, 2, 1, 2]].';
y = TimeDependentFiltering(noisy_signal, {tf(1), P_2}, sample_index_to_filter_index_mat);
% figure; spectrum(signal)
figure; spectrum(y)
% sound(signal, fs)
sound(y, fs)
%% Retrospective Howling Detection:
enable_false_alarm_detection = 1;
enable_gain_control_coping = 0;
is_plot_flag = 1;
[HowlingFreqTable, ~] = RetrospectiveHowlingDetectionFunc(Conf, y, signal, y, is_plot_flag,...
                                                          enable_false_alarm_detection, enable_gain_control_coping);
% disp(HowlingFreqTable(1:10,:))
%% Analysis - Finding First Howling Detection:
SpecificHowlingFreqTable = HowlingFreqTable(abs(HowlingFreqTable.Freq - f_pole) < 50, :);
figure; histogram(SpecificHowlingFreqTable.Avg_Grad, 100)
detected_avg_grad_mean = mean(SpecificHowlingFreqTable.Avg_Grad);
detected_avg_grad_std = std(SpecificHowlingFreqTable.Avg_Grad);

first_detection_delay_time_from_noise = SpecificHowlingFreqTable.Time(1) - 1.5;
first_detection_delay_time_from_full_buffer = SpecificHowlingFreqTable.Time(1) - 1.222;

disp({...
    ['detected_avg_grad_mean: ' num2str(detected_avg_grad_mean)];...
    ['detected_avg_grad_std: ' num2str(detected_avg_grad_std)];...
    ['first_detection_delay_time_from_noise: ' num2str(first_detection_delay_time_from_noise)];...
    ['first_detection_delay_time_from_full_buffer: ' num2str(first_detection_delay_time_from_full_buffer)]})

