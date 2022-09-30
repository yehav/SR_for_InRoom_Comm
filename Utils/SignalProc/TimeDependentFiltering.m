function [output_signal] = TimeDependentFiltering(input_signal, filters_cell, sample_index_to_filter_index_mat)
% Filter 'input_signal' with changing filters from 'filters_cell'.

% Theory:
% -------
% Y(z)/X(z) = (b0 + b1*z^-1 + b2*z^-2 + ... + bp*z^-p)/(1 + a1*z^-1 + a2*z^-2 + ... + aq*z^-q)
% => Y(z) = X(z)*(b0 + b1*z^-1 + b2*z^-2 + ... + bp*z^-p) - Y(z)*(a1*z^-1 + a2*z^-2 + ... + aq*z^-q)
% => y[n] = x[n]*b0 + x[n-1]*b1 + x[n-2]*b2 + ... + x[n-p]*bp - (y[n-1]*a1 + y[n-2]*a2 + ... + y[n-q]*aq)
% So, when a filter is changed, we should consider:
% * x[n-p:n-1]
% * y[n-q:n-1]

% Author: Yehav Alkaher.

%% For Debug:
if nargin == 0
    input_signal = ones(1,25);
    fs = 16e3;
    P_1 = tf([4], [1 1], 1/fs, 'variable','z^-1');
    % [y1, zf] = filter(P.Numerator{1}, P.Denominator{1}, ones(1, 10))
    P_2 = tf([1 2 3], [4 5 6], 1/fs, 'variable','z^-1');
    P_3 = tf([7 2 3], [4 5 6], 1/fs, 'variable','z^-1');
    P_4 = tf([1 1 1 1], [1], 1/fs, 'variable','z^-1');
    P_5 = tf([1 1 1], [1], 1/fs, 'variable','z^-1');
    P_6 = tf(2*[1 1 1], [1], 1/fs, 'variable','z^-1');
    filters_cell = {P_1, P_2, P_3, P_4, P_5, P_6};
    
    % sample_index_to_filter_index_mat = [1, 1; 15, 2];
    sample_index_to_filter_index_mat = [1, 4; 9, 4; 14, 4];
    sample_index_to_filter_index_mat = [1, 5; 9, 6; 14, 5];
end
%%
output_signal = [];

%% Check Inputs
% input type
assert(iscell(filters_cell) && ~isempty(filters_cell))
if isempty(sample_index_to_filter_index_mat)
    sample_index_to_filter_index_mat = [1, 1];
end

% 'sample_index_to_filter_index_mat(:, 1)'
assert( length( size(sample_index_to_filter_index_mat) ) == 2 )
assert( size(sample_index_to_filter_index_mat, 2) == 2 )
assert( issorted(sample_index_to_filter_index_mat(:, 1)) )

% both columns of 'sample_index_to_filter_index_mat'
assert( all(min(sample_index_to_filter_index_mat, [], 1) >= 1))
if sample_index_to_filter_index_mat(1, 1) ~= 1
    sample_index_to_filter_index_mat = [[1, 1]; sample_index_to_filter_index_mat];
end

% 'sample_index_to_filter_index_mat(:, 2)' & 'filters_cell'
assert( all(mod(sample_index_to_filter_index_mat(:, 2), 1) == 0) )
assert( max(sample_index_to_filter_index_mat(:, 2)) <= length(filters_cell))
%% Simplify Filters
for filter_idx = 1:length(filters_cell)
    P_itr = filters_cell{filter_idx};
    P_denom = P_itr.Denominator{1};
    P_itr.Numerator{1} = P_itr.Numerator{1}/P_denom(1);
    P_itr.Denominator{1} = P_itr.Denominator{1}/P_denom(1);
    filters_cell{filter_idx} = P_itr;
end
%% Filtering Initialization:
prev_filter_idx = 0;
zi = [];
% output_signal = [];
%% Apply Filters
for change_idx = 1:size(sample_index_to_filter_index_mat,1)
    curr_sample_index = sample_index_to_filter_index_mat(change_idx, 1);
    curr_filter_idx = sample_index_to_filter_index_mat(change_idx, 2);
    if size(sample_index_to_filter_index_mat, 1) >= change_idx + 1
        next_filters_change_sample_index = sample_index_to_filter_index_mat(change_idx + 1, 1);
    else
        next_filters_change_sample_index = length(input_signal)+1;
    end
    
    P_curr = filters_cell{curr_filter_idx};
    if prev_filter_idx == 0
        P_prev = tf(1);% tf([1], [1], 1/fs, 'variable','z^-1');
    else
        P_prev = filters_cell{prev_filter_idx};
    end
    
    % Parse Current Filter:
    P_curr_num = P_curr.Numerator{1};
    curr_p = find(P_curr_num ~= 0, 1,'last') - 1;
    P_curr_num = P_curr_num(1:(curr_p+1));
    
    P_curr_denom = P_curr.Denominator{1};
    curr_q = find(P_curr_denom ~= 0, 1,'last') - 1;
    P_curr_denom = P_curr_denom(1:(curr_q+1));
    
    % Take cropped input-signal with a samples buffer
    curr_input_signal = input_signal(curr_sample_index : (next_filters_change_sample_index - 1));
    % Prepare filter
    P_curr_num_padded = [P_curr_num, zeros(1, max(0, length(zi) - curr_p))];
    % Filter
    % Ex-1:
    %  [y1, zf] = filter([1 1 2 0 0], [1], ones(1, 10), [0 0 7 2])
    % Ex-2:
    %  [y1, zf] = filter([1 1 1], [1], ones(1, 15), [])
    %  [y1, zf] = filter([1 1 1], [1], ones(1, 15), zf)
    % Ex-3:
    %  [y1, zf] = filter([1], [1 1 1], ones(1, 14), [])
    %  [y1, zf] = filter([1], [1 1 1], ones(1, 14), zf)
    [curr_output_signal, zf] = filter(P_curr_num_padded, P_curr_denom, curr_input_signal, zi);
    % Get relevant output
    output_signal = [output_signal, curr_output_signal];
    % For next iteration:
    zi = zf(1:find(zf ~= 0, 1,'last')).';
    if isempty(zi)
        zi = [];
    end
    prev_filter_idx = curr_filter_idx;
end

% output_signal_est = filter(P_curr_num_padded, P_curr_denom, input_signal);
% assert( sum(abs(output_signal - output_signal_est)) == 0 )

if nargin == 0
    % sample_index_to_filter_index_mat = [1, 5; 9, 6; 14, 5];
    output_signal_est = ...
        [1, 2, 3*ones(1, 9 - 1 - 2),...
        4, 5, 6*ones(1, 14 - 9 - 2),...
        5, 4, 3*ones(1, length(input_signal) + 1 - 14 - 2)];
    assert( sum(abs(output_signal - output_signal_est)) == 0 )
end
end
