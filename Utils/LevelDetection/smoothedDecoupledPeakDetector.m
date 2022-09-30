function [y_l, y_r] = smoothedDecoupledPeakDetector(x_l, y_l_old, y_r_old, alpha_a, alpha_r)
% This function provides a digital implementation of the Level Detection stage,
% via a smoothed decoupled Peak Detector.

% The Level Detection stage is used to provide a smooth representation of
% the signal's level, and *may be applied at various places* in the side-chain.

% The gradual change of gain is due to the attack and release times.
% Attack and release times provide a degree of control over
% how quickly a compressor acts.
% * Attack  - Rise - Increasing the Gain-Reduction by
%           decreasing the gain of the compressor to the level 
%           determined by the ratio (once the signal overshoots the threshold).
% * Release - Fall - Decreasing the Gain-Reduction by
%           bringing the gain back up to the normal level
%           (once the signal has fallen below the threshold).

% The attack and the release times are usually introduced through a
% smoothing detector filter, e.g. exponential moving average (EMA):
% y[n] = alpha * y[n - 1] + (1 - alpha) * x[ n], where alpha is the filter coefficient.

% The time constant tau is defined as the time it takes for this system
% to reach 1-1/e of its final value.
% Thus, alpha = exp(-1/(tau*fs)).
% On the opposite, tau = -1/(fs*log(alpha)).
% For example, from 1 to 0, after tau seconds the value is:
% 1 + (0 - 1)*(1 - 1/e) = 1/e.

% The inputs to this stage are:
% 1) x_l - input signal level.
% 2) y_l_old - preceding output signal level.
% 3) y_r_old - preceding output signal level of the release stage.
% 4) alpha_a - the filter coefficient that corresponds to the attack time (rise-time).
% 5) alpha_r - the filter coefficient that corresponds to the release time (fall-time).

% The outputs of this stage are:
% 1) y_l - output signal level.
% 2) y_r - output signal level of the release stage.

% The decoupled peak detector produces a measured time constant of
% approximately alpha_A + alpha_R.

% Author: Yehav Alkaher.

y_r_ema = alpha_r * y_r_old + (1 - alpha_r) * x_l;
y_r = max(x_l, y_r_ema);

y_l = alpha_a * y_l_old + (1 - alpha_a) * y_r;

end
%%
% fs = 1e3;% [Hz]
% t = 0:1/fs:1;
% x_l = 1*(t <= 0.3) + 0.25*(t > 0.3);
% tau_a = 50e-3;% [sec]
% tau_r = 100e-3;% [sec]
% alpha_a = exp(-1/(tau_a*fs)); alpha_r = exp(-1/(tau_r*fs));
% y_l_old = 0; y_r_old = 0;
% figure; hold on;
% for k = 1:length(t)
%     [y_l, y_r] = smoothedDecoupledPeakDetector(x_l(k), y_l_old, y_r_old, alpha_a, alpha_r);
%     y_l_old = y_l;
%     y_r_old = y_r;
%     plot(t(k), x_l(k), 'b*', t(k), y_l, 'g*')
% end
% hold off;
% legend('input signal', 'output signal')
% title({['Output of the smoothed decoupled peak detector circuit. ',...
%     'Attack Time = 50 ms, release time = 100 ms.'];...
%     'Inset, a close-up of the decoupled smoothed peak detectors as it responds to a sudden level change at 300 ms.';...
%     'Pay attention that the measured release time contant is 150 ms (\tau_A + \tau_R).'})
%%
% fs = 1e3;% [Hz]
% t = 0:1/fs:2;
% change_time = 1;% [sec]
% x_l = 1*(t <= change_time) + 0.25*(t > change_time);
% tau_a = 50e-3;% [sec]
% tau_r = 100e-3;% [sec]
% alpha_a = exp(-1/(tau_a*fs)); alpha_r = exp(-1/(tau_r*fs));
% y_l_old = 0; y_r_old = 0;
% y_l_old_ = 0; y_r_old_ = 0;
% figure; hold on;
% for k = 1:length(t)
%     [y_l, y_r] = smoothedDecoupledPeakDetector(x_l(k), y_l_old, y_r_old, alpha_a, alpha_r);
%     y_l_old = y_l;
%     y_r_old = y_r;
%     plot(t(k), x_l(k), 'b*', t(k), y_l, 'g*')
%     
%     [y_l_, y_r_] = smoothedDecoupledPeakDetector(-x_l(k), -y_l_old_, -y_r_old_, alpha_a, alpha_r);
%     y_l_old_ = -y_l_;
%     y_r_old_ = -y_r_;
%     plot(t(k), -y_l_, 'r*')
% end
% % Rise:
% plot(xlim, 1/exp(1)*[1, 1], 'k:')
% plot(tau_a*[1, 1]          , [0, 1-1/exp(1)], 'g--')
% plot((tau_a + tau_r)*[1, 1], [0, 1-1/exp(1)], 'r--')
% % Fall:
% plot(xlim, 0.25+0.75*1/exp(1)*[1, 1], 'k--')
% plot((change_time + tau_a + tau_r)*[1, 1], [0, 0.25+0.75*1/exp(1)], 'g--')
% plot((change_time + tau_a)*[1, 1]        , [0, 0.25+0.75*1/exp(1)], 'r--')
% 
% hold off;
% grid minor
% legend(...
%     'input signal',...
%     'output signal - smoothed decoupled peak detector',...
%     'output signal - reversed smoothed decoupled peak detector')
% title({['Output of the smoothed decoupled peak detector circuit. ',...
%     'Attack Time = 50 ms, release time = 100 ms.'];...
%     'Inset, a close-up of the decoupled smoothed peak detectors as it responds to a sudden level change at 300 ms.';...
%     'Pay attention that the measured release time contant is 150 ms (\tau_A + \tau_R).'})
