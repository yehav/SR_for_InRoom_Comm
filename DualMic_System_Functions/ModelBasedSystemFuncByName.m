function [m1, m2, y, measurements, HowlingFreqTable] = ModelBasedSystemFuncByName(...
    method_flag, Conf, u1, u2, k, g1, g2, filtersLength, noiseStd, rng_seed)
% This function returns the Performance-results of the matching Speech-Reinforcement system.
% 
% Used in:
% - DemonstrationScript.m
% 
% Author: Yehav Alkaher.

if nargin < 11
    speechFrameLength = [];
    LMS_stepSize = [];
elseif nargin < 10
    rng_seed = [];
end

if isempty(rng_seed)
    rng(123456);
else
    rng(rng_seed);
end

switch method_flag
    case 'SimpleGain'
        flag_MSD = 0;
        [m1, m2, y, measurements, HowlingFreqList] = MSD_GainControl_ModelBasedSystemFunc(Conf, u1, u2, k, g1, g2, filtersLength, noiseStd, flag_MSD);
    case 'MSD_GainControl'
        flag_MSD = 1;
        [m1, m2, y, measurements, HowlingFreqList] = MSD_GainControl_ModelBasedSystemFunc(Conf, u1, u2, k, g1, g2, filtersLength, noiseStd, flag_MSD);
    otherwise
        error('Choose ''method_flag''')
end

HowlingFreqTable = array2table(HowlingFreqList);
if ~isempty(HowlingFreqTable)
    HowlingFreqTable.Properties.VariableNames={...
        'Time', 'Freq', 'Magnitude',...
        'Avg_Grad', 'Avg_Grad_STD', 'Avg_2nd_Grad', 'Second_Grad_RMS',...
        'Percent_above_Zero', 'Percent_above_STD','Peak_Dominance_Ratio', 'RMS_val', 'Moving_RMS_Median'...
        };
end

end

