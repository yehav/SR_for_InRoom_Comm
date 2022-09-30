function [alpha_vec_est , R_0_to_p, a_vec_est, estimation_error] = AutocorrelationLinearPrediction(signal_buffer,p,plot_flag)
% This function computes the Linear Prediction coefficients, using the
% 'levinson' function, and returns the estimated model coefficient vector
% of length p.

% Input:
% - signal_buffer - a vector of the last N samples of y(n) - numeric vector.
% - p - order of the Linear Prediction model - integer scalar bigger than 0.
%       p is the number of coefficients for an autoregressive model of order p.
% - plot_flag - 1 or 0.
% Output:
% - alpha_vec_est
% - R_0_to_p
% - a_vec_est
% - error - the error corresponds to the variance of the white noise error.
% uses: xcorr, toeplitz, levinson.
% used in: SymmetricItakuraSaituDistanceCalc.
% Example:
% y(n) = 0.1*y(n-1) - 0.8*y(n-2) + w(n)
% H(z) = Y(Z) / X(Z) = 1 / A(z)
% A(z) = 1 - a1*z^(-1) - ... - ap*z^(-p)
%      = alpha_vec.' * [z^0 , z^(-1), ... , z^(-p)].'
% alpha_vec = [1 , -a1 , ... , -ap] = [1 , -a_vec]
% 
% Author: Yehav Alkaher.
%% Checks:
if ((p <= 0) || (mod(p,1) ~= 0))
    error('p is not an integer scalar bigger than 0');
end

if (~isvector(signal_buffer))
    error('signal_buffer is not a vector.')
end
%% Linear Prediction of order p - via Autocorrelation Method:
[r,lags] = xcorr(signal_buffer,'biased');
r(lags<0) = [];

[alpha_vec_est,e] = levinson(r,p);
if (plot_flag)
    disp(['error_by_levinson = ' num2str(e)])
    disp(['alpha_est = ' num2str(alpha_vec_est)])
end

a_vec_est = -alpha_vec_est;
a_vec_est(1) = 1;
if (plot_flag)
    disp(['a_est = ' num2str(a_vec_est)])
end
% pinv(toeplitz(r(1:p)))*r(2:p+1)

R_0_to_p = toeplitz(r(1:p+1));
estimation_error = alpha_vec_est * R_0_to_p * alpha_vec_est.';
if (plot_flag)
    disp(['error = ' num2str(estimation_error)])
end
% toeplitz(r(1:p))*a_vec_est(2:p+1).'-r(2:p+1)

all_pole_gain = sqrt(r(1:p+1).' * alpha_vec_est.');
% if ~isnan(all_pole_gain)
%     assert(abs(all_pole_gain - sqrt(estimation_error)) < 10^(floor(log10(all_pole_gain)) - 10))
% end
end
%% For tests:
% clc; close all; clear
% 
% alpha_vec = [1 0.1 -0.8];
% v = 0.4;
% w = sqrt(v)*randn(15000,1);
% x = filter(1,alpha_vec,w);
% signal_buffer = x(1:5000);
% 
% p = numel(alpha_vec)-1;
% 
% plot_flag = 1;
% [alpha_vec_est , R_0_to_p, a_vec_est, estimation_error] = AutocorrelationLinearPrediction(signal_buffer,p,plot_flag)
