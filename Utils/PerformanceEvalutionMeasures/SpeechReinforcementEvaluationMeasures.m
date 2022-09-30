function [measurements] = SpeechReinforcementEvaluationMeasures(reference_signal, distorted_signal, fs)
% This function extracts performance-evaluation measures for
% speech-reinforcement algorithms.
% 
% Inputs:
% *) reference_signal
% *) distorted_signal
% *) fs - sampling frequency - [Hz]
% Output:
% *) measurmeants - struct
% 
% Author: Yehav Alkaher.
%% Check Input:
if any(size(reference_signal) ~= size(distorted_signal))
    error('size(reference_signal) ~= size(distorted_signal)')
end
%% Initialize struct:
measurements = [];

%% Ineligibility Measure:
Overall_reference_to_distorted_STOI = stoi( reference_signal, distorted_signal, fs );
measurements.STOI = Overall_reference_to_distorted_STOI;

%% Quality Measure:
Overall_reference_to_distorted_PESQ = pesq_mex_fast_vec( reference_signal, distorted_signal, fs );
measurements.PESQ = Overall_reference_to_distorted_PESQ;

%% Distortion Measures:
ItakuraSaitu_num_of_coefficients = round(fs/1e3) + 2;% LinearPredictiveOrder:
OverallSymmISDist = SymmetricItakuraSaituDistanceCalc(reference_signal, distorted_signal, ItakuraSaitu_num_of_coefficients);
measurements.OverallSymmISDist = OverallSymmISDist;

[symm_IS_dist_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, @SymmetricItakuraSaituDistanceCalc);
st_SymmISDist_average = mean(symm_IS_dist_result_vec);
st_SymmISDist_median = median(symm_IS_dist_result_vec);
measurements.ShortTermSymmISDistAvg = st_SymmISDist_average;
measurements.ShortTermSymmISDistMedian = st_SymmISDist_median;

[IS_dist_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, @ItakuraSaituDistanceCalc);
st_ISDist_average = mean(IS_dist_result_vec);
st_ISDist_median = median(IS_dist_result_vec);
measurements.ShortTermISDistAvg = st_ISDist_average;
measurements.ShortTermISDistMedian = st_ISDist_median;

[IS_measure_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, @ItakuraSaituMeasureCalc);
st_ISMeas_average = mean(IS_measure_result_vec);
st_ISMeas_median = median(IS_measure_result_vec);
measurements.ShortTermISMeasAvg = st_ISMeas_average;
measurements.ShortTermISMeasMedian = st_ISMeas_median;

[LLR_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, @LogLikelihoodRatioMeasureCalc);
st_LLRMeas_average = mean(LLR_result_vec);
st_LLRMeas_median = median(LLR_result_vec);
measurements.ShortTermLLRMeasAvg = st_LLRMeas_average;
measurements.ShortTermLLRMeasMedian = st_LLRMeas_median;

[EffectiveGain_result_vec] = LongTermSpectralDistanceMeasure(reference_signal, distorted_signal, fs, @EffectiveGainMeasureCalc);
st_EffectiveGain_median = median(EffectiveGain_result_vec);
measurements.ShortTermEffectiveGainMedian = st_EffectiveGain_median;
end

