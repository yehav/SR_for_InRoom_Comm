function [MSD, GC] = GenerateMsdGainControlStructs(fs, desired_gain,...
                                                   MSD_Conf, enable_false_alarm_detection, enable_gain_control_coping)
% This function generates the structs, used for the MSD-GC segment.
%
% Input:
% *) fs - sampling frequncy
% *) desired_gain
% *) Optional - MSD_Conf
%               + calculation_frame_length_time
%               + num_of_detection_frames
%               + num_of_frames
%               + detect_by_msd_alone_flag (Relevant in 'MsdBasedHowlingDetect')
% *) Optional - enable_false_alarm_detection
%                * Relevant in 'MsdBasedHowlingDetect'.
% *) Optional - enable_gain_control_coping
%                * Relevant in 'MsdGainControl' & 'MsdBasedHowlingDetect'.
%                * In 'MsdGainControl'        - Reduce 'num_of_frames' in History-Buffer after a howling detection.
%                * In 'MsdBasedHowlingDetect' - Use a softer 'ratio_above_mid_mag_upper_th' (in 'Underdamped_FA_Detect_Thresholds')
%                                                until enough time passes from last howling detection.
% 
% Author: Yehav Alkaher.
%% Check Input:
if (nargin < 3) || isempty(MSD_Conf)
    MSD_Conf = [];
    MSD_Conf.num_of_frames = 120;
    MSD_Conf.num_of_detection_frames = 6;
    MSD_Conf.calculation_frame_length_time = 30e-3; % about a speech frame length.
    MSD_Conf.detect_by_msd_alone_flag = 0;
end
MSD_Conf_field_names = {'num_of_frames', 'num_of_detection_frames', 'calculation_frame_length_time', 'detect_by_msd_alone_flag'};
assert(all(contains(fieldnames(MSD_Conf), MSD_Conf_field_names)))
assert(all(contains(MSD_Conf_field_names, fieldnames(MSD_Conf))))

if nargin < 4
    enable_false_alarm_detection = 1;
end

if nargin < 5
    enable_gain_control_coping = 1;
end

%% MSD:
MSD = [];
MSD.Const.fs = fs;
MSD.Const.num_of_frames = MSD_Conf.num_of_frames;
MSD.Const.num_of_detection_frames = MSD_Conf.num_of_detection_frames;
MSD.Const.msd_frame_length = 2^round(log(MSD_Conf.calculation_frame_length_time * fs)/log(2)); % [samples]
MSD.Const.msd_frame_shift = ceil(fs*10e-3);

MSD.Const.minimal_howling_detection_delay = ...
    (MSD.Const.msd_frame_length +...
    MSD.Const.msd_frame_shift * (MSD.Const.num_of_detection_frames - 1))/fs;
disp(['minimal_howling_detection_delay = ' num2str(MSD.Const.minimal_howling_detection_delay)])

MSD.Const.howling_history_delay = ...
    (MSD.Const.msd_frame_length +...
    MSD.Const.msd_frame_shift * (MSD.Const.num_of_frames - 1))/fs;
disp(['howling_history_delay = ' num2str(MSD.Const.howling_history_delay)])

MSD.msd_buffer = nan + zeros(1, MSD.Const.msd_frame_length);
MSD.magnitude_history_mat = nan + zeros(MSD.Const.msd_frame_length/2, MSD.Const.num_of_frames);

MSD.Detect_Method.detect_by_msd_alone_flag = MSD_Conf.detect_by_msd_alone_flag;

%% - Minimal howl-energy-threshold:
MSD.Const.safety_below_minimal_howl_energy_th = 5;% dB
% - Uniform Threshold:
% % MSD.Const.minimal_howl_energy_threshold = -95;% dB
% - Using "Equal-Loudness Contours" - choosing a 'phon' value:
N = MSD.Const.msd_frame_length/2;
freqs = 0:(0.5*fs/N):0.5*fs*(1-1/N);
minimal_phon = 20;
[spl] = iso226_2003_spl_from_loudness_lvl(minimal_phon,freqs,true);% dB SPL
reference_sound_pressure_air = 20e-6;
MSD.Const.minimal_howl_energy_threshold = spl(:) + db(reference_sound_pressure_air) - MSD.Const.safety_below_minimal_howl_energy_th;

MSD.Const.minimal_phon = minimal_phon;
MSD.Const.reference_sound_pressure_air = reference_sound_pressure_air;

% t = 0:(1/fs):(10000 - 1)/fs;
% mag_in_dB = MSD.Const.minimal_howl_energy_threshold;
% amplitude_val = 10^(mag_in_dB/20);
% m2 = amplitude_val*sin(2*pi*2e3*t);
% sound(m2, fs)

%% - Howl FA Detection
MSD.FA_Detect.Enable = enable_false_alarm_detection;% <-- Relevant in 'MsdBasedHowlingDetect'.
MSD.FA_Detect.Enable_GainControl_Coping = enable_gain_control_coping;
% ^ Implications:
%   * In 'MsdGainControl'        - Reduce 'num_of_frames' in History-Buffer after a howling detection.
%   * In 'MsdBasedHowlingDetect' - Use a softer 'ratio_above_mid_mag_upper_th' (in 'Underdamped_FA_Detect_Thresholds')
%                                   until enough time passes from last howling detection.
MSD.FA_Detect.min_num_of_frames_since_last_detection = 60;

MSD.FA_Detect.low_rms_threshold = 200;
MSD.FA_Detect.valid_rms_threshold = 500;
MSD.FA_Detect.sequence_size_threshold = 5;

%% --- Underdamped_Fine_FA_Detect_Thresholds:
Und_FA_Detect_Thresholds = [];
Und_FA_Detect_Thresholds.EndingLowRMS.Label_21_23_percent_above_median_of_avg_grad_th = 0.7;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_11_grad_lower_boundary = -20;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_1_grad_lower_boundary = -20;

Und_FA_Detect_Thresholds.EndingHighRMS.Label_22_percent_above_std_lower_th = 0.9;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_1_grad_lower_boundary = -20;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_24_percent_above_std_lower_th = 0.74;

Und_FA_Detect_Thresholds.momentary_gradient_estimation_boundary = 25;
Und_FA_Detect_Thresholds.ratio_above_mid_mag_upper_th = 0.1;

MSD.FA_Detect.Underdamped_Fine_FA_Detect_Thresholds = Und_FA_Detect_Thresholds;

%% --- Underdamped_MoreFine_FA_Detect_Thresholds:
Und_FA_Detect_Thresholds = [];
Und_FA_Detect_Thresholds.EndingLowRMS.Label_21_23_percent_above_median_of_avg_grad_th = 0.7;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_11_grad_lower_boundary = -20;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_1_grad_lower_boundary = -20;

Und_FA_Detect_Thresholds.EndingHighRMS.Label_22_percent_above_std_lower_th = 0.95;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_1_grad_lower_boundary = -20;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_24_percent_above_std_lower_th = 0.85;

Und_FA_Detect_Thresholds.momentary_gradient_estimation_boundary = 25;
Und_FA_Detect_Thresholds.ratio_above_mid_mag_upper_th = 0.15;
Und_FA_Detect_Thresholds.AfterDetection.ratio_above_mid_mag_upper_th = 0.2;

MSD.FA_Detect.Underdamped_MoreFine_FA_Detect_Thresholds = Und_FA_Detect_Thresholds;

%% --- Underdamped_Strict_FA_Detect_Thresholds:
Und_FA_Detect_Thresholds = [];
Und_FA_Detect_Thresholds.EndingLowRMS.Label_21_23_percent_above_median_of_avg_grad_th = 0.7;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_11_grad_lower_boundary = -5;
Und_FA_Detect_Thresholds.EndingLowRMS.Label_1_grad_lower_boundary = -5;

Und_FA_Detect_Thresholds.EndingHighRMS.Label_22_percent_above_std_lower_th = 0.95;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_1_grad_lower_boundary = 0;
Und_FA_Detect_Thresholds.EndingHighRMS.Label_24_percent_above_std_lower_th = 0.86;

Und_FA_Detect_Thresholds.momentary_gradient_estimation_boundary = 25;
Und_FA_Detect_Thresholds.ratio_above_mid_mag_upper_th = 0.05;

MSD.FA_Detect.Underdamped_Strict_FA_Detect_Thresholds = Und_FA_Detect_Thresholds;

%% --- Increasing_FA_Detect_Thresholds:
Inc_FA_Detect_Thresholds = [];
Inc_FA_Detect_Thresholds.percent_all_above_zero_th = 0.5;

Inc_FA_Detect_Thresholds.EndingLowRMS.Label_400_percent_above_std_lower_th = 0.95;
Inc_FA_Detect_Thresholds.EndingLowRMS.Label_401_percent_above_median_of_avg_grad_th = 0.9;
Inc_FA_Detect_Thresholds.EndingLowRMS.Label_404_405_grad_boundary = -5;

Inc_FA_Detect_Thresholds.EndingHighRMS.Label_395_grad_boundary = 10;

MSD.FA_Detect.Increasing_FA_Detect_Thresholds = Inc_FA_Detect_Thresholds;

%% Level Detection:
GC = [];
%- for reversed smoothed-decoupled peak-detector:
gain_tau_a = 5e-3;% Fall-Time [sec]
gain_tau_r = 10e-3;% Rise-Time [sec]
% -> Effective Rise-Time is: gain_tau_a + gain_tau_r

GC.LevelDetect.Const.gain_alpha_a = exp(-1/(gain_tau_a * fs));
GC.LevelDetect.Const.gain_alpha_r = exp(-1/(gain_tau_r * fs));

GC.LevelDetect.gain_r_old = 0;
GC.LevelDetect.gain_curr = 1;
GC.LevelDetect.gain_dest = desired_gain;
GC.LevelDetect.gain_max = desired_gain;
GC.LevelDetect.available_desired_gain = nan;
GC.LevelDetect.gain_max_lower_boundary = 0.5;

round_to_upper_half = @(x) max(round(x), round(x + 0.5) - 0.5);
GC.LevelDetect.level_detection_start_time = 0*round_to_upper_half(MSD.Const.howling_history_delay)*fs;

GC.Time.Const.gain_control_renew_time = 1*fs;
GC.Time.Const.gain_control_release_time = 100e-3*fs;
GC.Time.Const.howling_detect_release_time = 60e-3*fs;
GC.Time.n_howling = nan;% time of last howling detection
GC.Time.n_next_renew = nan;% time of last 'renew_time'
end

