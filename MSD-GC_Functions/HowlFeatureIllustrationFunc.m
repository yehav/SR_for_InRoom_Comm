function [] = HowlFeatureIllustrationFunc(features, frequency_magnitude_threshold_dB, minimized_plot_flag, momentary_howl_avg_grad)
% Illustrate the Features of a history-buffer at a specific frequency-bin.
% 
% Author: Yehav Alkaher.

if nargin < 3
    minimized_plot_flag = 0;
    momentary_howl_avg_grad = [];
end
if nargin < 4
    momentary_howl_avg_grad = [];
end

moving_magnitude_avg_vec = features.MovingFilters.moving_magnitude_avg_vec;
sequences_of_low_rms_vec = features.LowRMS.sequences_of_low_rms_vec;
sequences_of_high_but_valid_rms_vec = features.LowRMS.sequences_of_high_but_valid_rms_vec;
% sequences_of_valid_rms_vec = features.LowRMS.sequences_of_valid_rms_vec;
num_of_low_rms_sequences = features.LowRMS.num_of_low_rms_sequences;
low_rms_start_indices = features.LowRMS.low_rms_start_indices;
median_of_avg_grad_per_low_rms_vec = features.LowRMS.MovingFilters.median_of_avg_grad_per_low_rms_vec;

if minimized_plot_flag
    figure('Position', [633.67 41.667 460.67 599.33]);
else
    figure;
end

p1 = plot(moving_magnitude_avg_vec);
hold on;
moving_magnitude_avg_vec_ = moving_magnitude_avg_vec(~isinf(moving_magnitude_avg_vec));
% p2 = plot(min(moving_magnitude_avg_vec_) + (max(moving_magnitude_avg_vec_) - min(moving_magnitude_avg_vec_))*sequences_of_valid_rms_vec, 'r');
p3 = plot(min(moving_magnitude_avg_vec_) + (max(moving_magnitude_avg_vec_) - min(moving_magnitude_avg_vec_))*sequences_of_high_but_valid_rms_vec, 'ro');
p4 = plot(min(moving_magnitude_avg_vec_) + (max(moving_magnitude_avg_vec_) - min(moving_magnitude_avg_vec_))*sequences_of_low_rms_vec, 'g');
for itr_idx = 1:num_of_low_rms_sequences
    median_of_avg_grad_per_low_rms_sequence = [num2str(median_of_avg_grad_per_low_rms_vec(itr_idx), 5) ' [dB/sec]'];
    text(low_rms_start_indices(itr_idx)+1, max(moving_magnitude_avg_vec_)+0.01*range(moving_magnitude_avg_vec_), median_of_avg_grad_per_low_rms_sequence)
end
p5 = plot(xlim, frequency_magnitude_threshold_dB*[1 1], 'k--');
hold off
y_limits = [min(moving_magnitude_avg_vec_), max(moving_magnitude_avg_vec_)] +...
            0.03*range(moving_magnitude_avg_vec_)*[-1, 1];
y_limits = [min(y_limits(1), frequency_magnitude_threshold_dB - 10),...
            max(y_limits(2), frequency_magnitude_threshold_dB + 10)];
ylim(y_limits)

h_legend = legend([p1, p4, p3, p5],'Moving Magnitude-Avg vec', 'Sequences of Low-RMS vec', 'Sequences of High-but-Valid-RMS vec', 'Frequency Magnitude-Threshold [dB]');
h_legend.Position(1) = 0.2;
h_legend.Position(2) = 0.6;

title_text_cell = {'Median of ''Average 1st-Gradient'' values in Low-RMS Sequence';...
                   'illustrated over the Magnitude Moving-Average'};
if ~isempty(momentary_howl_avg_grad)
    title_text_cell{3} = ['Howling Momentary Avg. 1st-Grad.: ' num2str(momentary_howl_avg_grad) ' [dB/sec]'];
end
title(title_text_cell)
end

