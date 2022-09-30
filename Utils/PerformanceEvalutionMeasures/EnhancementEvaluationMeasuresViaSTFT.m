function [matched_effective_gain, normalized_distorion] = EnhancementEvaluationMeasuresViaSTFT(input_signal, enhanced_signal, fs, desired_gain_factor, gain_resolution)
% This function extracts performance-evaluation measures that estimate the
% effective applied 'gain-factor' and the resulting 'normalized-distorion',
% using the STFT of the input and enhanced signals.
% 
% Input:
% *) input_signal
% *) enhanced_signal
% *) fs
% *) desired_gain_factor
% *) gain_resolution
% Output:
% *) matched_effective_gain - the effective applied-gain over time.
% *) normalized_distorion - the normalized signal-distortion.
% 
% Author: Yehav Alkaher.
%% Check Input:
if nargin < 5
    gain_resolution = 0.05;
end
%% Initialization:
gain_grid = 0 : gain_resolution : desired_gain_factor;
calculation_frame_length_time = 30e-3; % about twice a speech frame length.
frame_length = 2^round(log(calculation_frame_length_time * fs)/log(2)); % [samples]

if any(size(input_signal) ~= size(enhanced_signal))
    error('size(input_signal) ~= size(enhanced_signal)')
end

[input_s, input_f, input_t] = spectrogram(input_signal, frame_length, frame_length/2, [], fs,'yaxis');
[enhanced_s, enhanced_f, enhanced_t] = spectrogram(enhanced_signal, frame_length, frame_length/2, [], fs,'yaxis');

if any(input_f ~= enhanced_f) || any(input_t ~= enhanced_t)
    error('any(input_f ~= enhanced_f) || any(input_t ~= enhanced_t)')
end
%%
absolute_diff_values_per_gain = abs(repmat(enhanced_s(:),1,length(gain_grid)) - input_s(:)*gain_grid);
[least_maximal_absolute_diff_value, gain_grid_index] = min(max(absolute_diff_values_per_gain,[],1));
matched_effective_gain = gain_grid(gain_grid_index);
normalized_distorion = least_maximal_absolute_diff_value / max(abs(enhanced_s(:)));
end

