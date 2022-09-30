function [matched_effective_gain, normalized_distorion] = EnhancementEvaluationMeasuresViaMSE(input_signal, enhanced_signal, fs, desired_gain_factor, gain_resolution)
% This function extracts performance-evaluation measures that estimate the
% effective applied 'gain-factor' and the resulting 'normalized-distorion',
% using the MSE of the input and enhanced signals.
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

if any(size(input_signal) ~= size(enhanced_signal))
    error('size(input_signal) ~= size(enhanced_signal)')
end

%%
MSE_per_gain = sum((repmat(enhanced_signal(:),1,length(gain_grid)) - input_signal(:)*gain_grid).^2,1);
[least_MSE, gain_grid_index] = min(max(MSE_per_gain,[],1));
matched_effective_gain = gain_grid(gain_grid_index);
normalized_distorion = least_MSE / norm(enhanced_signal)^2;
end

