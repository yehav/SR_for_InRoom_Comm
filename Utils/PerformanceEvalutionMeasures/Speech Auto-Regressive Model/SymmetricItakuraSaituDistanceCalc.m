function [Symmetric_IS_dist] = SymmetricItakuraSaituDistanceCalc(signal_buffer_a, signal_buffer_b, p)
% This function computes the Symmetric Itakura-Saito distance between two signal vectors.

% The Symmetric Itakura-Saito distance is a measure of
% the difference between an original spectrum and an approximation of that
% spectrum.

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

if (~isvector(signal_buffer_a) || ~isvector(signal_buffer_b))
    error('the two input signal buffers must be vectors of the same size.')
end

% if (size(signal_buffer_a) ~= size(signal_buffer_b))
%     error('the two input signal buffers must be vectors of the same size.')
% end

%% Linear Prediction of order p - via Autocorrelation Method:
plot_flag = 0;
[alpha_vec_est , R_0_to_p_a, ~, ~] = AutocorrelationLinearPrediction(signal_buffer_a,p,plot_flag);
[beta_vec_est , R_0_to_p_b, ~, ~] = AutocorrelationLinearPrediction(signal_buffer_b,p,plot_flag);

%% Symmetric Itakura-Saitu Distance:
alpha_error = alpha_vec_est * R_0_to_p_a * alpha_vec_est.';
difference_error_a = (alpha_vec_est - beta_vec_est) * R_0_to_p_a * (alpha_vec_est - beta_vec_est).';

IS_dist_a = difference_error_a/alpha_error;


beta_error = beta_vec_est * R_0_to_p_b * beta_vec_est.';
difference_error_b = (alpha_vec_est - beta_vec_est) * R_0_to_p_b * (alpha_vec_est - beta_vec_est).';

IS_dist_b = difference_error_b/beta_error;


Symmetric_IS_dist = 0.5 * (IS_dist_a + IS_dist_b);


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
% [Symmetric_IS_dist] = SymmetricItakuraSaituDistanceCalc(signal_buffer_a, signal_buffer_b, p)
