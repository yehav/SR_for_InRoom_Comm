% This script implements the initialization for the
% Speech Reinforcement System's model-based functions.
% 
% Author: Yehav Alkaher.

u1 = u1(:);
%% Initializations:
fs = Conf.fs;

fftPlot = Conf.Functions.fftPlot;
spectrum = Conf.Functions.spectrum;

EnableNoise = Conf.EnableNoise;

isPlotFlag = Conf.isPlotFlag;
% isPlotFlag = 1;% adaptive filters' impulse responses, final spectrums.
% isPlotFlag = 2;% final spectrums, Error plot, impulse response comparison, ERLE.
% isPlotFlag = 3;% h_ESF_NRF effect in the fourier domain.
isDispFlag = Conf.isDispFlag;

if (EnableNoise == 1)
    disp('Enabling Noise');
end
%% Delays:
% SystemDelayBlock - Models the processing delay of the system, since
%              receiving the signal from the microphone to emitting it from
%              the loudspeaker.
SystemTimeDelay = Conf.Delays.systemTimeDelay;% [msec]
SystemDelayBlock = ceil(SystemTimeDelay*fs);%[samples]
disp(['the System DelayBlock is: ' num2str(SystemDelayBlock) ' [samples]']);
% MinPropDelayBlock - Models the minimal electro-acoustic delay of the loudspeaker to
%              microphone path, we would like to answer.
%              This delay time is lower bounded by the adaptive filter multiplication processing delay.
MinPropagationTimeDelay = Conf.Delays.minPropagationTimeDelay;% [msec]
MinPropDelayBlock = ceil(MinPropagationTimeDelay*fs);%[samples]
disp(['the Min Propagation DelayBlock is: ' num2str(MinPropDelayBlock) ' [samples]']);

%% Parameters and Vectors:
% 1. Result vectors:
m1 =  zeros(length(u1),1);
m2 =  zeros(length(u1),1);
f1 = zeros(length(u1),1);
f2 = zeros(length(u1),1);
y =  zeros(length(u1),1);


% 2. Assisting buffers:
bufferLength = max(SystemDelayBlock + MinPropDelayBlock + filtersLength , SystemDelayBlock + length(g1));
y_buffer = zeros(1,bufferLength);

% For direct response, we take the last 'length(g1)' samples of the buffer:
RirIndices = (bufferLength-SystemDelayBlock-length(g1)+1) : (bufferLength-SystemDelayBlock) ;

% 3. Howling Detection:
HowlingFreqList = [];

