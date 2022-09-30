% This script implements the initialization for
%   Speech Reinforcement System scripts.
% 
% Author: Yehav Alkaher.


% Properties:
c = 340;% Sound velocity [m/s]

% Sampling Frequency:
fs = 16e3;%22.05e3;%44.1e3;

% Random Seed:
noise_rng_seed = 123456;

% Microphone Saturation - Max Amplitude before Clipping:
max_signal_amplitude = 1;

% HPF to remove DC from loudspeaker's output:
f_cut = 20;% [Hz]
remove_DC_filter = fir1(2^13,f_cut/(fs/2),'high');
remove_DC_filter = remove_DC_filter/max(remove_DC_filter);
% freqz(remove_DC_filter,1,length(remove_DC_filter))

% FFT plot properties:
fftPlot = @(sig) plot(linspace(-1,1,length(sig)),abs(fftshift(fft(sig))));

% 'mySpectrogram' Properties:
timeFrameLength = 0.03;% [sec]
frameLength_estimated = timeFrameLength * fs;
speechFrameLength = 2^round(log(frameLength_estimated)/log(2));%128; % [samples]
overlap_ratio = 0.5; % [%]
newFigureFlag = 0; % 1 for new figure.
spectrum = @(sig) mySpectrogram(sig, speechFrameLength, overlap_ratio, fs, newFigureFlag);

% STFT properties:
% Conf.STFT.N = 256;% Num of frequency bins.
% Conf.STFT.L = 256;% Num of timeframe samples.
% Conf.STFT.D = 64;% Num of timeshift samples.

%  Configuration:
Conf.fs = fs;
Conf.max_signal_amplitude = max_signal_amplitude;
Conf.remove_DC_filter = remove_DC_filter;
Conf.Functions.fftPlot = fftPlot;
Conf.Functions.spectrum = spectrum;
Conf.LinearPredictiveOrder = floor(fs/1e3)+2;