function [m1, m2, y, measurements, HowlingFreqList] = MSD_GainControl_ModelBasedSystemFunc(Conf, u1, u2, k, g1, g2, filtersLength, noiseStd, flag_MSD)
% This script implements the
%   Dual-Microphone Speech Reinforcement System with Howling-Control
%    for In-Room Speech Communication.

% Inputs:
% u1 - near-end speech - source input signal in the microphone m1.
% u2 - near-end speech - source input signal in the microphone m2.
% g1 - RIR from loudspeaker to micrphone m1 - Loudspeaker-Enclosure-Microphone (LEM) path.
% g2 - RIR from loudspeaker to micrphone m2 - Loudspeaker-Enclosure-Microphone (LEM) path.
% w - white noise from loudspeaker.
% b - background noise to microphone m1.
% K - gain factor by the user.

% Additional Buffers:
% m1 - speaker's microphone signal.
% m1 - howling detection microphone signal.
% y - far-end speech - loudspeaker output signal.
% f1 - feedback echo signal from loudspeaker to microphone m1.
% f2 - feedback echo signal from loudspeaker to microphone m2.

% Author: Yehav Alkaher.

tic

Initialize_model_based_functions
spectrum = Conf.Functions.spectrum;

%% Gain Control:
% User Configured Gain-Control:
% Check if 'k' changes along time:
configured_gain = k;
maximal_configured_gain = max(k);
changing_gain_flag = 0;
if length(k) == length(u1)
    changing_gain_flag = 1;
    configured_gain = k(1);
else
    assert(length(k) == 1)
end

[MSD, GC] = GenerateMsdGainControlStructs(fs, configured_gain);

Tracked_values.gain_values = [];

measurements = [];
%% System Core:
for n = 1:length(u1)
    %% Input Vectors:
    m1(n) = MicrophoneSignalAcquisition(u1(n), f1(n), noiseStd, EnableNoise, Conf.max_signal_amplitude);
    m2(n) = MicrophoneSignalAcquisition(u2(n), f2(n), noiseStd, EnableNoise, Conf.max_signal_amplitude);
    
    MSD.msd_buffer = [MSD.msd_buffer(2:end), m2(n)];
    
    %% Gain to be Applied:
    if changing_gain_flag
        if k(n) ~= configured_gain
            GC.Time.n_next_renew = n-1;
            GC.Time.n_howling = n-1 - GC.Time.Const.gain_control_renew_time;
        end
        configured_gain = k(n);
    end
    %% Gain-Control:
    if (flag_MSD == 1)
        [MSD, GC, HowlingFreqList] = MsdGainControl(n, fs, configured_gain, MSD, GC, Conf, HowlingFreqList, isDispFlag, isPlotFlag, m1);
        
        Tracked_values.gain_values = [Tracked_values.gain_values, GC.LevelDetect.gain_curr];
    end
    
    %% Output:
    w_n = EnableNoise*noiseStd*randn(1);
    
    if (flag_MSD == 1)
        y_n = GC.LevelDetect.gain_curr * m1(n) + w_n;
    else
        y_n = configured_gain*m1(n) + w_n;
    end
    
    y_buffer = [y_buffer(2:end) , y_n];
    
    f1_n = FIRFilterResponse( y_buffer(RirIndices) , g1 );
    f2_n = FIRFilterResponse( y_buffer(RirIndices) , g2 );
    
    %% Results Check:
    if (isnan(y_n) || isnan(f1_n) || isnan(f2_n))
        error('Nan values exist!')
    end
    if (isinf(y_n) || isinf(f1_n) || isinf(f2_n))
        error('Inf values exist!')
    end
    %% Update Result Vectors:
    y(n) = y_n;
    if (n < length(u1))
        f1(n+1) = f1_n;
        f2(n+1) = f2_n;
    end
    
    %% End of Calculations.
end
%% System Core - End.
%% Retrospective Measurements:
if (isPlotFlag ~= 0)
    figure('Position', get(0, 'Screensize'));
    subplot(2,2,1); spectrum(u1); title('$u_1$', 'Interpreter','latex');
    
    subplot(2,2,2);
    u1_fft_magnitude = abs(fftshift(fft(u1)));
    plot(linspace(-fs/2, fs/2, length(u1_fft_magnitude)), u1_fft_magnitude);
    title('FFT$\{ u_1 \}$', 'Interpreter','latex');
    ylabel('Magnitude')
    xlabel('Frequency [Hz]')
    
    subplot(2,2,3); spectrum(y); title('$y$', 'Interpreter','latex');
    
    subplot(2,2,4);
    y_fft_magnitude = abs(fftshift(fft(y)));
    plot(linspace(-fs/2, fs/2, length(y_fft_magnitude)), y_fft_magnitude);
    title('FFT$\{ y \}$', 'Interpreter','latex');
    ylabel('Magnitude')
    xlabel('Frequency [Hz]')
    
    
    if (flag_MSD == 1)
        figure('Position', get(0, 'Screensize'));
        subplot(4,1,1); spectrum(u1); title('u1');
        subplot(4,1,2); spectrum(y); title('y');
        subplot(4,1,3); plot(0:1/fs:(length(u1)-1)/fs, Tracked_values.gain_values);
        xlim([0, (length(y)-1)/fs]); ylim([0 maximal_configured_gain * 1.1]); colorbar;
        title('gain values over time:');
        ylabel('Gain')
        xlabel('Time [sec]')
        
        subplot(4,1,4);
        if ~isempty(HowlingFreqList)
            plot(HowlingFreqList(:,1), HowlingFreqList(:,2), '*');
        end
        xlim([0, (length(y)-1)/fs]); ylim([0, fs/2]); colorbar;
        title('Howling frequencies over time:');
        ylabel('Howling Frequency [Hz]')
        xlabel('Time [sec]')
        
        
        figure;
        ax1 = subplot(2,1,1); spectrum(y);
        title('y');
        hold on
        if ~isempty(HowlingFreqList)
            plot(HowlingFreqList(:,1), HowlingFreqList(:,2)/1e3, 'ro', 'MarkerSize', 10);
        end
        colorbar('off')
        xlabel(''); xticklabels('')
        xlim([0, (length(y)-1)/fs]);
        
        ax2 = subplot(2,1,2); plot(0:1/fs:(length(u1)-1)/fs, Tracked_values.gain_values);
        title('gain values over time:');
        xlim([0, (length(y)-1)/fs]); ylim([0 maximal_configured_gain * 1.1]);
        grid on
        ylabel('Gain')
        xlabel('Time [sec]')
        
        p1 = get(ax1, 'Position');
        p2 = get(ax2, 'Position');
        p1(2) = p2(2)+p2(4);
        set(ax1, 'pos', p1);
        
        linkaxes([ax1 ax2],'xy')
    end
end
%%

%% Return Measurements:
[measurements] = SpeechReinforcementEvaluationMeasures(u1, y, fs);
if (flag_MSD == 1)
    measurements.RMS_applied_gain = rms_calc(Tracked_values.gain_values.', length(u1), 1);
end
[stft_matched_effective_gain, stft_normalized_distorion] = EnhancementEvaluationMeasuresViaSTFT(u1, y, fs, maximal_configured_gain);
measurements.STFT_matched_effective_gain = stft_matched_effective_gain;
measurements.STFT_normalized_distorion = stft_normalized_distorion;
[mse_matched_effective_gain, mse_normalized_distorion] = EnhancementEvaluationMeasuresViaMSE(u1, y, fs, maximal_configured_gain);
measurements.MSE_matched_effective_gain = mse_matched_effective_gain;
measurements.MSE_normalized_distorion = mse_normalized_distorion;

toc
end

%% Inner Functions:
function [m_n] = MicrophoneSignalAcquisition(u_n, f1_n, micNoiseStd, EnableNoise, max_input_amplitude)
% This function is responsible of acquiring the microphone signals, which
% include: Source signal, Background Noise, and Echo originated from the loudspeaker.
% Inputs:
%   1) u1_n - Source signal from the microphone at time 'n' - scalar.
%           - near-end speech - source input signal in the microphone m1.
%   2) f1_n - Feedback echo signal from loudspeaker to the microphone at time 'n' - scalar.
%   4) micNoiseStd - Standard Deviation of microphone thermal noise - scalar.
%   5) EnableNoise - Enable noise flag - 0 or 1.
% Outputs:
%   1) m_n - microphone signal at time 'n' - scalar.
b_n = micNoiseStd*randn(1);
m_n = EnableNoise*b_n + u_n + f1_n;
if abs(m_n) > max_input_amplitude
    m_n = sign(m_n) * max_input_amplitude;
    warning('abs(m_n) > max_input_amplitude')
end
end