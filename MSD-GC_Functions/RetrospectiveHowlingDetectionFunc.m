function [HowlingFreqTable, MSD] = RetrospectiveHowlingDetectionFunc(Conf, m2, u1, y, is_plot_flag,...
                                                                     enable_false_alarm_detection, enable_gain_control_coping, GC_coping_bypass_flag,...
                                                                     MSD_Conf)
% This script performs a Retrospective Howling Detection Post-Factum.
% 
% For Debug:
% 1)
% t = 0:(1/fs):(3*1024 - 1)/fs; m2 = 10*sin(2*pi*2e3*t);
% max(abs(fft(m2)))/(length(t)/2)
% [psd_buffer, zero2piRange] = PSD_calc(m2, m2, 0, 0, 0);
% psd_buffer_normalized = psd_buffer/((length(psd_buffer)/2)^2);
% mag_db = 10*log10(max(psd_buffer)) - 20*log10(length(psd_buffer)/2);
% disp([num2str(10^(mag_db/20)) ' - ' num2str(10) ' == 0 : ' num2str(abs(10^(mag_db/20) - 10) < 1e-6)])
% 2)
% [psd_buffer, zero2piRange] = PSD_calc(m2, m2, 0, 0, 0);
% disp(['20*log10(10) == ' num2str(10*log10(max(psd_buffer)) - 20*log10(length(psd_buffer)/2))]);
% figure; spectrogram(m2, 1024, 512, [], fs,'yaxis');
% 3)
% RetrospectiveHowlingDetectionFunc(Conf, m2, 0, 0, 0)
% 
% Author: Yehav Alkaher.

tic

spectrum = Conf.Functions.spectrum;
fs = Conf.fs;

%% Check Input:
if nargin < 6
    enable_false_alarm_detection = 1;
end

if nargin < 7
    enable_gain_control_coping = 1;
end

if nargin < 8
    GC_coping_bypass_flag = 0;
end

if nargin < 9
    MSD_Conf = [];
end

%% Gain Control:
% MSD:
[MSD, GC] = GenerateMsdGainControlStructs(fs, 0, MSD_Conf, enable_false_alarm_detection, enable_gain_control_coping);
% Optional - Comment for deminished detection:
GC.Time.Const.howling_detect_release_time = 0;
% Optional - Bypass-Mode for gain-control coping conditions:
if GC_coping_bypass_flag
    MSD.Const.num_of_frames = MSD.FA_Detect.min_num_of_frames_since_last_detection;
    MSD.magnitude_history_mat = nan + zeros(MSD.Const.msd_frame_length/2, MSD.Const.num_of_frames);
    % Condition in 'MsdBasedHowlingDetect' applies:
    % Und_Thresholds.ratio_above_mid_mag_upper_th = Und_Thresholds.AfterDetction.ratio_above_mid_mag_upper_th;
end

HowlingFreqList = [];

%% System Core:
for n = 1:length(m2)
    %% Input Vectors:
    MSD.msd_buffer = [MSD.msd_buffer(2:end), m2(n)];
    
    %% Gain-Control:
    disp_details_flag = 0;
    [MSD, GC, HowlingFreqList] = MsdGainControl(n, fs, 0, MSD, GC, Conf, HowlingFreqList, disp_details_flag, 0, []);
    
    %% End of Calculations.
end
%% System Core - End.

HowlingFreqTable = array2table(HowlingFreqList);
if ~isempty(HowlingFreqTable)
    HowlingFreqTable.Properties.VariableNames={...
        'Time', 'Freq', 'Magnitude',...
        'Avg_Grad', 'Avg_Grad_STD', 'Avg_2nd_Grad', 'Second_Grad_RMS',...
        'Percent_above_Zero', 'Percent_above_STD','Peak_Dominance_Ratio', 'RMS_val', 'Moving_RMS_Median'...
        };
end
%% Retrospective Analysis:
if ((is_plot_flag ~= 0) && ~isempty(HowlingFreqTable))
    figure;
    ax1 = subplot(2,1,1); spectrum(u1);
    hold on
    if contains(ax1.XLabel.String, 'mins')
        plot(HowlingFreqTable.Time/60, HowlingFreqTable.Freq/1e3, 'ro', 'MarkerSize', 10);
    else
        plot(HowlingFreqTable.Time, HowlingFreqTable.Freq/1e3, 'ro', 'MarkerSize', 10);
    end
    colorbar('off')
    title('(a)');
    title('$u_1(n)$ --- Acquired Speaker''s Speech', 'Interpreter','latex');
    ax2 = subplot(2,1,2); spectrum(y);
    hold on
    if contains(ax2.XLabel.String, 'mins')
        plot(HowlingFreqTable.Time/60, HowlingFreqTable.Freq/1e3, 'ro', 'MarkerSize', 10);
    else
        plot(HowlingFreqTable.Time, HowlingFreqTable.Freq/1e3, 'ro', 'MarkerSize', 10);
    end
    colorbar('off')
    title('$y(n)$ --- Loudspeaker''s Output', 'Interpreter','latex')
    linkaxes([ax1 ax2],'xy')
    
    % [HowlingFreqTable, MSD_retro] = RetrospectiveHowlingDetectionFunc(Conf, m2, u1, y, 1)
    % saveas(gcf,[targetFolder 'spectrogram_comparison_with_howling_detecion.eps'],'epsc')
    % saveas(gcf,[targetFolder 'spectrogram_comparison_with_howling_detecion.png'],'png')
    % savefig(gcf,[targetFolder 'spectrogram_comparison_with_howling_detecion.fig'])
    % close;
end
%%
% sound(y,fs)
toc
end