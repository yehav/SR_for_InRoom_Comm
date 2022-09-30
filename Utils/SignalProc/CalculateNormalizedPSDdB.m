function normalized_psd_zero2pi_db = CalculateNormalizedPSDdB(buffer)
% This function calculates the common periodogram (power spectral density
% estimate) of x and y.
% Input - x and y are vectors.
%
% Output - psd_xy, zero2piRange.
%       For normalized result:
%       *) psd_buffer_normalized = psd_buffer/((length(psd_buffer)/2)^2);
%       *) mag_db = 10*log10(max(psd_buffer_normalized));
%
% Test:
% fs = 16e3; t = 0:(1/fs):(1024 - 1)/fs;
% freq = 128*fs/length(t);
% mag_in_dB = 60;
% amplitude_val = 10^(mag_in_dB/20);
% if mod(length(t),fs/freq) ~= 0
%     error('frequency bin is not exact')
% end
% m2 = amplitude_val*sin(2*pi*freq*t);
% normalized_psd_zero2pi_db = CalculateNormalizedPSDdB(m2);
%
% mag_db = max(normalized_psd_zero2pi_db);
% disp(['PSD magnitude per Hz [dB]: ' num2str(mag_db) ' ~ ' num2str(mag_in_dB)]);
%
% figure; plot(normalized_psd_zero2pi_db);
%
% Author: Yehav Alkaher.

psd_buffer = abs(fft(buffer(:))).^2;
% The fft indices moves from [ [0 to pi) , [-pi to 0) ] or [0,2*pi).
% Take only [0,pi] or [0,pi).
nfft = length(psd_buffer);
zero2piRange = 1:ceil(nfft/2);%floor(nfft/2)+1;

psd_buffer_normalized = psd_buffer/((length(psd_buffer)/2)^2);
normalized_psd_zero2pi_db = 10*log10(psd_buffer_normalized(zero2piRange));

end

