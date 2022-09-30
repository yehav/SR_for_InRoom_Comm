function [spectral_measure_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, spectral_measure_function)
% This function calculated the Short-Term
% Spectral Distance Measures Based on LPC, such as:
% * Symmetric Itakura-Saito (IS) distance measure
% * Itakura-Saito (IS) measure
% * log-likelihood ratio (LLR) measure
% for time-frames of the input reference- and distorted-signals.
% 
% Input:
% *) reference_signal
% *) distorted_signal
% *) fs - sampling frequency - [Hz]
% *) spectral_measure_function - a function that receives the arguments: (reference_frame, distorted_frame, order)
% Output:
% *) spectral_measure_result_vec
% 
% Author: Yehav Alkaher.
%%
frame_length_for_measurments = round(fs * 30e-3);%[Hz*sec = samples]
% order - LinearPredictiveOrder - ItakuraSaitu_num_of_coefficients:
order = round(fs/1e3) + 2;

reference_buffers = buffer(reference_signal, frame_length_for_measurments, floor(frame_length_for_measurments/2), 'nodelay');
distorted_buffers = buffer(distorted_signal, frame_length_for_measurments, floor(frame_length_for_measurments/2), 'nodelay');

if 0 ~= mod(length(reference_signal) - ceil(frame_length_for_measurments/2), frame_length_for_measurments)/frame_length_for_measurments
    reference_buffers = reference_buffers(:,1:end-1);
end
distorted_buffers = distorted_buffers(:,1:size(reference_buffers,2));

spectral_measure_result_vec = cellfun(@(reference_frame, distorted_frame) spectral_measure_function(reference_frame, distorted_frame, order),...
    mat2cell(reference_buffers, size(reference_buffers,1), ones(size(reference_buffers,2),1)),...
    mat2cell(distorted_buffers, size(distorted_buffers,1), ones(size(distorted_buffers,2),1)));

% Assert that nan values come from zero-arrays:
nan_value_indices = find(isnan(spectral_measure_result_vec));
assert(...
    all(all(reference_buffers(:, nan_value_indices) == 0)) ||...
    all(all(distorted_buffers(:, nan_value_indices) == 0)))

spectral_measure_result_vec = spectral_measure_result_vec(~isnan(spectral_measure_result_vec));
end

