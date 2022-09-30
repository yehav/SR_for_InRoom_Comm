function [x_rms] = rms_calc(X, num_of_samples, num_of_channels)
% This function provides a simple RMS calculation on blocks of samples.
% This function work for a mono and stereo samples.
% For stereo - a mono down-mix (L+R)/2 of the programme audio is used.

% The input:
% 1) X - matrix of size num_of_samples x num_of_channels.
% 2) num_of_samples - number of samples in the frame.
% 3) num_of_channels - number of stero channels.

% The output: x_rms - Root-Mean-Sqaure value of the input-signal's mono down-mix.

% Author: Yehav Alkaher.

if any(size(X) ~= [num_of_samples, num_of_channels])
    error('size(X) ~= [num_of_samples, num_of_channels]')
end

x_mono_vec = mean(X, 2);
x_rms = sqrt(mean(x_mono_vec.^2));
end

