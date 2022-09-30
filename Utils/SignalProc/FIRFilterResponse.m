function [output] = FIRFilterResponse(x_buffer,h)
% This function calculates the output a 'causal' impulse response on a
% buffer.

% Inputs:
% x_buffer = [x(n-p+1), x(n-p+2), ... , x(n)] - row vector of length p.
% h = [h(0), h(1), ... , h(p)] - row vector of length p.

% Output:
% output =
% [x(n), x(n-1) , ... , x(n-p+1)] * transpose( [h(0), h(1), ... , h(p)] )

% Author: Yehav Alkaher.

if( size(x_buffer) ~= size(h) )
    disp('size(x_buffer) ~= size(h)')
    return;
end

x_buffer = x_buffer(:).';

output = fliplr(x_buffer) * h.';
end

