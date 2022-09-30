function [filtered_logical_vec] = DetectSequenceOfIndications(logical_vec, sequence_size_threshold, plot_flag, pad_intro_flag)
% Detect a sequence of 'True' logical-values within a logical vector.

% Input:
% - logical_vec
% - sequence_size_threshold - minimal number of 'True' logical-values in row
% - plot_flag
% - pad_intro_flag

% Output:
% - filtered_logical_vec - new 'logical_vec' without the filtered indications.

% For Debug:
% logical_vec = [1   1   1   1   1   1   1   1   1   1   0   0   0   0,...
%     0   0   1   0   0   0   0   0   0   0   0   0   0   1,...
%     1   1   0   0   0   0   0   0   0   0   0   0   0   0,...
%     0   0   0   1   1   1   1   1   1   1   1   1];
% sequence_size_threshold = 3;
% plot_flag = 1;
% pad_intro_flag = 0;
% 
% [filtered_logical_vec] = DetectSequenceOfIndications(logical_vec, sequence_size_threshold, plot_flag, pad_intro_flag)
%
% 
% Author: Yehav Alkaher.

if nargin < 3
    plot_flag = 0;
    pad_intro_flag = 0;
end
if nargin < 4
    pad_intro_flag = 0;
end

assert(sequence_size_threshold > 1)

work_logical_vec = [pad_intro_flag*ones(1, sequence_size_threshold-1), logical_vec(:).'];
moving_average_vec = conv(work_logical_vec, ones(1, sequence_size_threshold)/sequence_size_threshold, 'valid');

sequence_indices = find(moving_average_vec == 1);

fixed_sequence_indices = sequence_indices(:) + (-1*(sequence_size_threshold-1):0);
fixed_sequence_indices = fixed_sequence_indices(fixed_sequence_indices(:) > 0);
fixed_sequence_indices = unique(fixed_sequence_indices);

filtered_logical_vec = 0*logical_vec;
filtered_logical_vec(fixed_sequence_indices) = 1;
filtered_logical_vec = (filtered_logical_vec > 0);

if plot_flag == 1
    figure;
    plot(logical_vec, 'k')
    hold on;
    stem(filtered_logical_vec, 'b')
    hold off;
end

end

