function [] = mySpectrogram(signal, frameLength, overlap_ratio, fs, newFigureFlag)
% % properties:
% frameLength = 128; % [samples]
% overLap = 0.5; % [%]
% fs = 1e3; % [Hz]
% newFigureFlag is 1 for new figure.

% soundsc(signal,fs);

if newFigureFlag == 1
    figure;
end

spectrogram(signal, frameLength, round(overlap_ratio*frameLength), [], fs,'yaxis');

end

