function [EffectiveGain_measure] = EffectiveGainMeasureCalc(reference_signal, distorted_signal, p)
% This function computes the Effective-Gain between two signal vectors.

% The Effective-Gain is measured as the ratio between the all-pole gains
% of the clean and enhanced signals.

% Input:
% - reference_signal - a vector of the last N samples of y_a(n) - numeric vector.
% - distorted_signal - a vector of the last N samples of y_b(n) - numeric vector.
% - p - order of the Linear Prediction model - integer scalar bigger than 0.
%       p is the number of coefficients for an autoregressive model of order p.
% Output:
% - error - the error corresponds to the variance of the white noise error.
% uses: AutocorrelationLinearPrediction.
% used in: - used as a distortion measure of speech signals.
% 
% Author: Yehav Alkaher.
%% Checks:
if ((p <= 0) || (mod(p,1) ~= 0))
    error('p is not an integer scalar bigger than 0');
end

if (~isvector(reference_signal) || ~isvector(distorted_signal))
    error('the two input signal buffers must be vectors of the same size.')
end

% if (size(signal_buffer_a) ~= size(signal_buffer_b))
%     error('the two input signal buffers must be vectors of the same size.')
% end

%% Linear Prediction of order p - via Autocorrelation Method:
plot_flag = 0;
[~, ~, ~, all_pole_gain_squared_reference] = AutocorrelationLinearPrediction(reference_signal,p,plot_flag);
[~, ~, ~, all_pole_gain_squared_distorted] = AutocorrelationLinearPrediction(distorted_signal,p,plot_flag);

%% Effective-Gain Measure:
EffectiveGain_measure = sqrt(all_pole_gain_squared_distorted / all_pole_gain_squared_reference);

end
%% For tests:
% clc; close all; clear
% 
% alpha_vec = [1 0.1 -0.8];
% v = 0.4;
% w = sqrt(v)*randn(15000,1);
% x = filter(1,alpha_vec,w);
% signal_buffer_a = x(1:5000);
% signal_buffer_b = x(501:5500);
% 
% p = numel(alpha_vec)-1;
% 
% plot_flag = 1;
% [IS_measure] = ItakuraSaituMeasureCalc(signal_buffer_a, signal_buffer_b, p)
% [Symmetric_IS_dist] = SymmetricItakuraSaituDistanceCalc(signal_buffer_a, signal_buffer_b, p)
% [LLR_measure] = LogLikelihoodRatioMeasureCalc(signal_buffer_a, signal_buffer_b, p)