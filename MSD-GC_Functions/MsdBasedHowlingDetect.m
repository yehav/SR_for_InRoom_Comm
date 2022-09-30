function [nBin, candidates, measurements] = MsdBasedHowlingDetect(MSD, fs)
% msdEvaluate Calculates the Magnitude Slope Deviation of all data points
%in magnitude data buffer. Uses double-differentiation to
%identify straight-line magnitude growth (dB scale) of component
%frequencies.

% Input:
% - MSD struct - contains:
%     * magnitude_history_mat
%     * msd_frame_length
%     * msd_frame_shift
%     * num_of_frames
% - fs - sampling frequency [Hz]

% Output:
% - nBin - all possible frequency bins (if peak-picking not active, algorithms search all bins).
% - candidates - indices of suspected frequency-howl bins.
% - measurements - of candidate frequency-howls.
% 
% Author: Yehav Alkaher.
%%
assert(size(MSD.magnitude_history_mat,1) == MSD.Const.msd_frame_length/2)
assert(MSD.FA_Detect.Enable_GainControl_Coping || size(MSD.magnitude_history_mat,2) == MSD.Const.num_of_frames)

num_of_frames = size(MSD.magnitude_history_mat, 2);
assert(num_of_frames >= MSD.Const.num_of_detection_frames)

detection_frame_indices = (num_of_frames - MSD.Const.num_of_detection_frames + 1):num_of_frames;
magnitude_history_mat = MSD.magnitude_history_mat(:, detection_frame_indices);
msd_frame_length = MSD.Const.msd_frame_length;
msd_frame_shift = MSD.Const.msd_frame_shift;
%%
% i) Peak-Picking:
nBin = 1:msd_frame_length/2;% if peak-picking not active, algorithms search all bins
% grad = diff(MSD.magnitude_history_mat,1,2) / (msd_frame_shift / fs);% first derivative
% nBin = nBin(mean(grad ,2) < 0); % dB
magnitude_history_mat = magnitude_history_mat(nBin, :);
%%
% Time-Diff between Frames:
dt = msd_frame_shift / fs;

% ii) Check MSD:
grad = diff(magnitude_history_mat,1,2) / dt;% first derivative - dB/sec
grad_mean = mean(grad, 2);
grad_std = std(grad, [], 2);

% grad_change - diff(grad,1,2) / (msd_frame_shift / fs)
grad_change = diff(magnitude_history_mat,2,2) / dt^2;% second derivative - dB/sec^2
grad_change_mean = mean(grad_change, 2);
% grad_change_std = std(grad_change, [], 2);

% For an unbiased estimator (assuming 2nd derivative is 0 for a frequency-howl),
% the RMSD is the square root of the variance, known as the standard deviation.
grad_change_rms = sqrt(mean(abs(grad_change).^2,2)); % dB/sec^2

%%

% figure;
% scatter(grad_mean, abs(grad_change_mean), [], (nBin-1)*fs/msd_frame_length, 'filled', 'MarkerFaceAlpha', 0.3);
% cb = colorbar;
% cb.Label.String = 'Howling Frequency [Hz]';
% xlabel('Magnitude Slope [dB/sec]')
% ylabel('Magnitude Slope Deviation [dB/sec^2]')
%%
% iii) Detect Howls:
if MSD.Detect_Method.detect_by_msd_alone_flag
    candidate_indices = find((grad_change_rms < 3.5e3));% & (grad_change_rms > 50)); % dB/sec^2
else
    % iii.A) Increasing Howls:
    increasing_bin_indices = find(grad_mean > 0);
    candidate_increasing_howl_indices = ...
        (grad_std(increasing_bin_indices) < 200) &...
        (abs(grad_change_mean(increasing_bin_indices))  < 2e3) &...
        (grad_change_rms(increasing_bin_indices) < 1e4);
    
    % iii.B) Underdamped Howls:
    underdamped_bin_indices = find(grad_mean < 0);
    candidate_underdamped_howl_indices = ...
        (grad_mean(underdamped_bin_indices) > -1000) &...
        (grad_std(underdamped_bin_indices) < 200) &...
        (abs(grad_change_mean(underdamped_bin_indices)) < 1e4) &...
        (grad_change_rms(underdamped_bin_indices) < 3e4);
    
    % iv) Combine howling candidates (among nBin):
    candidate_indices = [...
        increasing_bin_indices(candidate_increasing_howl_indices);...
        underdamped_bin_indices(candidate_underdamped_howl_indices)];
end

candidates = nBin(candidate_indices);

if isempty(candidates)
    measurements = [];
    return
end

%% Measurement Reports:
% (candidates - 1)*fs/msd_frame_length
mean_energey_per_candidate = mean(magnitude_history_mat(candidate_indices,:), 2);

howl_history_vectors_mat = MSD.magnitude_history_mat(candidates,:);
howl_history_1st_grad_mat = diff(howl_history_vectors_mat,1,2) / dt;
centered_howl_history_1st_grad_mat = howl_history_1st_grad_mat - min(0, grad_mean(candidate_indices));

% parameters:
percent_above_zero = mean(centered_howl_history_1st_grad_mat >=0, 2);
percent_above_std = mean(centered_howl_history_1st_grad_mat + grad_std(candidate_indices) >=0, 2);
peak_dominance_ratio = max(centered_howl_history_1st_grad_mat,[],2) ./ ( max(centered_howl_history_1st_grad_mat,[],2) - min(centered_howl_history_1st_grad_mat,[],2));
rms_val = sqrt(mean(abs(centered_howl_history_1st_grad_mat).^2,2));

grad_window_length = MSD.Const.num_of_detection_frames - 1;
moving_rms_mat = sqrt(conv2(abs(centered_howl_history_1st_grad_mat).^2, ones(1, grad_window_length)/grad_window_length, 'valid'));
moving_rms_median_vec = median(moving_rms_mat,2);
% moving_rms_std_vec = std(moving_rms_mat, [], 2);
moving_grad_avg_mat = conv2(howl_history_1st_grad_mat, ones(1, grad_window_length)/grad_window_length, 'valid');
assert(all(moving_grad_avg_mat(:, end) - grad_mean(candidate_indices) < 1e-7))


measurements = [...
    (candidates.'-1)*fs/msd_frame_length,...
    mean_energey_per_candidate,...
    ...
    grad_mean(candidate_indices),...
    grad_std(candidate_indices),...
    abs(grad_change_mean(candidate_indices)),...
    grad_change_rms(candidate_indices),...
    ...
    percent_above_zero,...
    percent_above_std,...
    peak_dominance_ratio,...
    rms_val,...
    moving_rms_median_vec,...
    ...
    []];

measurementsTable=array2table(measurements);
measurementsTable.Properties.VariableNames={...
    'Freq', 'Magnitude',...
    'Avg_Grad', 'Avg_Grad_STD', 'Avg_2nd_Grad', 'Second_Grad_RMS',...
    'Percent_above_Zero', 'Percent_above_STD','Peak_Dominance_Ratio', 'RMS_val', 'Moving_RMS_Median'...
    };
% disp(measurementsTable)

%% False-Alarm Detection
howl_label_vec = [];
if MSD.FA_Detect.Enable
    % Und_Thresholds = MSD.FA_Detect.Underdamped_Fine_FA_Detect_Thresholds;
    % Und_Thresholds = MSD.FA_Detect.Underdamped_Strict_FA_Detect_Thresholds;
    Inc_Thresholds = MSD.FA_Detect.Increasing_FA_Detect_Thresholds;
    
    % - iso226_2003_spl_from_loudness_lvl:
    % minimal_phon = MSD.Const.minimal_phon;
    % reference_sound_pressure_air = MSD.Const.reference_sound_pressure_air;
    % freqs = (candidates.'-1)*fs/msd_frame_length;
    % frequency_magnitude_threshold_dB_vec = iso226_2003_spl_from_loudness_lvl(minimal_phon,freqs,true) + db(reference_sound_pressure_air);
    % assert( all(MSD.Const.minimal_howl_energy_threshold(candidates) == frequency_magnitude_threshold_dB_vec) )
    
    howl_label_vec = zeros(length(candidate_indices), 1);
    for candidate_idx = 1:length(candidate_indices)
        % assert(candidates(candidate_idx) == nBin(candidate_indices(candidate_idx)))
        freq = (candidates(candidate_idx)-1)*fs/msd_frame_length;
        howl_history_vec = howl_history_vectors_mat(candidate_idx, :);
        howl_avg_grad = grad_mean(candidate_indices(candidate_idx));
        howl_avg_grad_std = grad_std(candidate_indices(candidate_idx));
        frequency_magnitude_threshold_dB = MSD.Const.minimal_howl_energy_threshold(candidates(candidate_idx));
        
        % Feature Extraction:
        [features] = HowlFeatureExtraction(MSD, howl_history_vec, howl_avg_grad, howl_avg_grad_std, fs);
        if 0
            minimized_plot_flag = 1;
            HowlFeatureIllustrationFunc(features, frequency_magnitude_threshold_dB, minimized_plot_flag, howl_avg_grad);
        end
        
        % False-Alarm Detection:
        if freq > 1500% Hz
            % Und_Thresholds = MSD.FA_Detect.Underdamped_Fine_FA_Detect_Thresholds;
            Und_Thresholds = MSD.FA_Detect.Underdamped_MoreFine_FA_Detect_Thresholds;
            
            if MSD.FA_Detect.Enable_GainControl_Coping
                if num_of_frames <= MSD.FA_Detect.min_num_of_frames_since_last_detection
                    % Use a softer 'ratio_above_mid_mag_upper_th'
                    %   until enough time passes from last howling detection.
                    Und_Thresholds.ratio_above_mid_mag_upper_th = Und_Thresholds.AfterDetection.ratio_above_mid_mag_upper_th;
                end
            end
            
            Und_Thresholds = rmfield(Und_Thresholds, 'AfterDetection');
        else
            Und_Thresholds = MSD.FA_Detect.Underdamped_Strict_FA_Detect_Thresholds;
        end
        
        show_details_flag = 0;
        [Label, ~] = HowlingEvaluation_mex(freq, howl_avg_grad, features, frequency_magnitude_threshold_dB, Und_Thresholds, Inc_Thresholds, show_details_flag);
        
        howl_label_vec(candidate_idx) = Label;
    end
    
    howl_label_vec = howl_label_vec .* ((candidates.'-1)*fs/msd_frame_length > 1e3);
    
    candidates = candidates(howl_label_vec > 0 );
    measurements = measurements(howl_label_vec > 0, :);
end
%%
end

