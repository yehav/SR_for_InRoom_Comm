function [features] = HowlFeatureExtraction(MSD, howl_history_vec, howl_avg_grad, howl_avg_grad_std, fs)
% This script extracts features for further howling False-Alarm detection.
% 
% Author: Yehav Alkaher.
%% Properties From MSD:
assert(fs == MSD.Const.fs)
msd_frame_shift = MSD.Const.msd_frame_shift;
dt = msd_frame_shift / fs;

% Post-Detection Configuration:
low_rms_threshold = MSD.FA_Detect.low_rms_threshold;
valid_rms_threshold = MSD.FA_Detect.valid_rms_threshold;
sequence_size_threshold = MSD.FA_Detect.sequence_size_threshold;

%% main vectors:
howl_history_vec_1st_grad = diff(howl_history_vec,1,2) / dt;
howl_history_vec_2nd_grad = diff(howl_history_vec,2,2) / dt^2;
centered_howl_history_vec_1st_grad = howl_history_vec_1st_grad - howl_avg_grad;% min(0, howl_avg_grad);

features = [];

features.howl_history_vec_1st_grad = howl_history_vec_1st_grad;
% features.howl_history_vec_2nd_grad = howl_history_vec_2nd_grad;
% features.centered_howl_history_vec_1st_grad = centered_howl_history_vec_1st_grad;

%% parameters:
features.percent_all_above_zero = mean(howl_history_vec_1st_grad >=0);

% Making sure there is no lower decay-rate in the sequence:
% features.percent_above_zero = mean(centered_howl_history_vec_1st_grad >=0);
features.percent_above_std = mean(centered_howl_history_vec_1st_grad + howl_avg_grad_std >=0);

% More Parameters:
% features.peak_dominance = max(centered_howl_history_vec_1st_grad) / ( max(centered_howl_history_vec_1st_grad) - min(centered_howl_history_vec_1st_grad));
% features.rms_val = sqrt(mean(abs(centered_howl_history_vec_1st_grad).^2,2));

% Moving Filters:
mag_window_length = MSD.Const.num_of_detection_frames;
grad_window_length = mag_window_length - 1;
% features.MovingFilters.mag_window_length = mag_window_length;
% features.MovingFilters.grad_window_length = grad_window_length;
% - Moving Centered-Gradient RMS:
moving_rms_vec = sqrt(conv(abs(centered_howl_history_vec_1st_grad).^2, ones(1, grad_window_length)/grad_window_length, 'valid'));
features.MovingFilters.moving_centered_grad_rms_vec = moving_rms_vec;
% features.MovingFilters.moving_centered_grad_rms_median = median(moving_rms_vec);
% features.MovingFilters.moving_centered_grad_rms_std = std(moving_rms_vec);
% - Moving Gradient Average:
features.MovingFilters.moving_grad_avg_vec = conv(howl_history_vec_1st_grad, ones(1, grad_window_length)/grad_window_length, 'valid');
% assert(features.MovingFilters.moving_grad_avg_vec(end) - howl_avg_grad < 1e-7)
% - Moving Magnitude Average:
% moving_magnitude_avg_vec = conv(howl_history_vec(2:end), ones(1, grad_window_length)/grad_window_length, 'valid');
features.MovingFilters.moving_magnitude_avg_vec = conv(howl_history_vec, ones(1, mag_window_length)/mag_window_length, 'valid');
%% extra parameters - squences_of_low_rms_vec & their 1st gradient:
sequences_of_low_rms_vec = DetectSequenceOfIndications(moving_rms_vec < low_rms_threshold, sequence_size_threshold);
sequences_of_valid_rms_vec = DetectSequenceOfIndications(moving_rms_vec < valid_rms_threshold, sequence_size_threshold);
sequences_of_high_but_valid_rms_vec = (~sequences_of_low_rms_vec) & sequences_of_valid_rms_vec;
% Divide to Sub-Sequences:
low_rms_start_indices = find(diff(sequences_of_low_rms_vec) == 1) + 1;
low_rms_end_indices = find(diff(sequences_of_low_rms_vec) == -1);
is_ending_with_low_rms = 0;
if sequences_of_low_rms_vec(1) == 1
    low_rms_start_indices = [1, low_rms_start_indices];
end
if sequences_of_low_rms_vec(end) == 1
    low_rms_end_indices = [low_rms_end_indices, length(sequences_of_low_rms_vec)];
    is_ending_with_low_rms = 1;
end
assert(length(low_rms_start_indices) == length(low_rms_end_indices))

num_of_low_rms_sequences = length(low_rms_start_indices);

% - Sanity-Check:
sequences_sum = 0;
for itr_idx = 1:num_of_low_rms_sequences
    sequences_sum = sequences_sum + sum(sequences_of_low_rms_vec(low_rms_start_indices(itr_idx):low_rms_end_indices(itr_idx)));
end
assert(sequences_sum == sum(sequences_of_low_rms_vec))

% Moving Filters:
moving_grad_avg_vec = features.MovingFilters.moving_grad_avg_vec;
moving_magnitude_avg_vec = features.MovingFilters.moving_magnitude_avg_vec;

mean_of_avg_grad_per_low_rms_vec = zeros(1, num_of_low_rms_sequences) + nan;
median_of_avg_grad_per_low_rms_vec = zeros(1, num_of_low_rms_sequences) + nan;
median_of_magnitude_avg_per_low_rms_vec = zeros(1, num_of_low_rms_sequences) + nan;
for itr_idx = 1:num_of_low_rms_sequences
    mean_of_avg_grad_per_low_rms_vec(itr_idx) = mean(moving_grad_avg_vec(low_rms_start_indices(itr_idx):low_rms_end_indices(itr_idx)));
    median_of_avg_grad_per_low_rms_vec(itr_idx) = median(moving_grad_avg_vec(low_rms_start_indices(itr_idx):low_rms_end_indices(itr_idx)));
    median_of_magnitude_avg_per_low_rms_vec(itr_idx) = median(moving_magnitude_avg_vec(low_rms_start_indices(itr_idx):low_rms_end_indices(itr_idx)));
end

% Conclude:
features.LowRMS.sequences_of_low_rms_vec = sequences_of_low_rms_vec;
% features.LowRMS.sequences_of_valid_rms_vec = sequences_of_valid_rms_vec;
features.LowRMS.sequences_of_high_but_valid_rms_vec = sequences_of_high_but_valid_rms_vec;

features.LowRMS.low_rms_start_indices = low_rms_start_indices;
features.LowRMS.low_rms_end_indices = low_rms_end_indices;
features.LowRMS.is_ending_with_low_rms = is_ending_with_low_rms;
features.LowRMS.num_of_low_rms_sequences = num_of_low_rms_sequences;

% features.LowRMS.MovingFilters.mean_of_avg_grad_per_low_rms_vec = mean_of_avg_grad_per_low_rms_vec;
features.LowRMS.MovingFilters.median_of_avg_grad_per_low_rms_vec = median_of_avg_grad_per_low_rms_vec;
% features.LowRMS.MovingFilters.median_of_magnitude_avg_per_low_rms_vec = median_of_magnitude_avg_per_low_rms_vec;
