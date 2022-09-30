function [IS_measure] = ItakuraSaituMeasureCalc(reference_signal, distorted_signal, p)
% This function computes the Itakura-Saito distance between two signal vectors.

% The Itakura-Saito measure (or Itakura-Saito divergence) is a measure of
% the difference between an original spectrum and an approximation of that
% spectrum.
% Although it is not a perceptual measure, it is intended to reflect
% perceptual (dis)similarity.

% Input:
% - signal_buffer_a - a vector of the last N samples of y_a(n) - numeric vector.
% - signal_buffer_b - a vector of the last N samples of y_b(n) - numeric vector.
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
[alpha_vec_est , R_0_to_p_a, ~, all_pole_gain_squared_a] = AutocorrelationLinearPrediction(reference_signal,p,plot_flag);
[beta_vec_est , ~, ~, all_pole_gain_squared_b] = AutocorrelationLinearPrediction(distorted_signal,p,plot_flag);

%% Itakura-Saitu Measure:
all_pole_gain_a = sqrt(all_pole_gain_squared_a);
all_pole_gain_b = sqrt(all_pole_gain_squared_b);

IS_measure = (all_pole_gain_a/all_pole_gain_b)*((beta_vec_est * R_0_to_p_a * beta_vec_est.')/(alpha_vec_est * R_0_to_p_a * alpha_vec_est.')) + log(all_pole_gain_b/all_pole_gain_a) - 1;

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
