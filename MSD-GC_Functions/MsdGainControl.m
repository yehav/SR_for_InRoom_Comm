function [MSD, GC, HowlingFreqList] = MsdGainControl(n, fs, desired_gain, MSD, GC, Conf, HowlingFreqList, disp_details_flag, is_plot_flag, m1)
% This function implements the Magnitude-Slope-Deviation (MSD) based
% Gain-Control Segment.
% 
% Input:
% *) n - sample index - integer
%
% *) fs - sampling frequency - double
% 
% *) desired_gain
% 
% *) MSD - struct:
%          a) Const - struct:
%                     a.1) num_of_frames    - integer
%                     a.2) msd_frame_length - integer
%                     a.3) msd_frame_nfft   - integer
%                     a.4) msd_frame_shift  - integer
%          b) msd_buffer - array of length: msd_frame_length
%          c) magnitude_history_mat - array of size: msd_frame_nfft x num_of_frames
% 
% *) GC  - struct:
%          a) LevelDetect - struct:
%                    a.1) Const - struct:
%                          a.1.i)   gain_alpha_a - integer
%                          a.1.ii)  gain_alpha_r - integer
%                    a.2) gain_r_old
%                    a.3) gain_curr
%                    a.4) gain_dest
%                    a.5) gain_max
%          b) Time - struct:
%                    b.1) Const - struct:
%                          b.1.i)   gain_control_renew_time     - integer
%                          b.1.ii)  gain_control_release_time   - integer
%                          b.1.iii) howling_detect_release_time - integer
%                    b.2) n_howling - integer
% 
% *) Conf - configuration struct
% 
% *) HowlingFreqList - array
% 
% *) disp_details_flag - 0 (don't display) or 1 (display)
% 
% *) is_plot_flag - 0, 1 (don't plot) or 2 (plot)
% 
% *) m1 - in case 'is_plot_flag' is 1.
% 
% Output:
% *) MSD - struct
% *) GC  - struct
% *) HowlingFreqList - array
% 
% Author: Yehav Alkaher.

%% Check Input:
if nargin < 7
    HowlingFreqList = [];
    disp_details_flag = 0;
    is_plot_flag = 0;
end

spectrum = Conf.Functions.spectrum;
%% Set Constants:
msd_frame_length = MSD.Const.msd_frame_length;
msd_frame_shift = MSD.Const.msd_frame_shift;
num_of_frames = MSD.Const.num_of_frames;

minimal_howl_energy_threshold = MSD.Const.minimal_howl_energy_threshold;
is_single_minimal_howl_energy_threshold = (length(minimal_howl_energy_threshold(:)) == 1);
is_per_freq_minimal_howl_energy_threshold = (length(minimal_howl_energy_threshold(:)) == msd_frame_length/2);
assert(is_single_minimal_howl_energy_threshold | is_per_freq_minimal_howl_energy_threshold)

gain_alpha_a = GC.LevelDetect.Const.gain_alpha_a;
gain_alpha_r = GC.LevelDetect.Const.gain_alpha_r;

gain_control_release_time = GC.Time.Const.gain_control_release_time;
gain_control_renew_time = GC.Time.Const.gain_control_renew_time;
howling_detect_release_time = GC.Time.Const.howling_detect_release_time;

%% Gain-Control:
%% 1) Check For 'Release-Time' - to update 'gain_dest':
% Set 'desired_gain_available' - for gradual gain increment:
if n > GC.Time.n_next_renew
    GC.LevelDetect.available_desired_gain = min(desired_gain, GC.LevelDetect.gain_max + 1);
    
    if GC.LevelDetect.available_desired_gain < desired_gain
        GC.Time.n_next_renew = GC.Time.n_next_renew + gain_control_renew_time;
    end
end
% Set 'gain_max':
if (n > (GC.Time.n_howling + gain_control_renew_time))
    GC.LevelDetect.gain_max = GC.LevelDetect.available_desired_gain;
end
% Set 'gain_dest':
if (n > (GC.Time.n_howling + gain_control_release_time))
    GC.LevelDetect.gain_dest = GC.LevelDetect.gain_max;
end
%% 2) Check For Howling - to start reducing the gain:
if (n >= msd_frame_length) && ( mod( n - msd_frame_length , msd_frame_shift) == 0 )
    normalized_psd_zero2pi_db = CalculateNormalizedPSDdB(MSD.msd_buffer + db2mag(-600)*randn(size(MSD.msd_buffer)));
    zero2piRange = 1:length(normalized_psd_zero2pi_db);
    % figure; plot(normalized_psd_zero2pi_db)
    
    MSD.magnitude_history_mat = [MSD.magnitude_history_mat(:,2:end), normalized_psd_zero2pi_db]; % add new data
    
    if n >= (msd_frame_length + msd_frame_shift * (num_of_frames - 1))
        MSD_copy = MSD;
        
        flag_gain_control_coping_fail = 0;
        if MSD.FA_Detect.Enable_GainControl_Coping && ~isnan(GC.Time.n_next_renew)
            %- curr_howling_time
            curr_howling_time = n/fs;% seconds
            %- last_howling_time
            %last_howling_time = GC.Time.n_howling/fs;% seconds
            %- last_gain_update_time
            if n > GC.Time.n_next_renew
                last_gain_update_time = GC.Time.n_next_renew/fs;
            else
                last_gain_update_time = (GC.Time.n_next_renew - gain_control_renew_time)/fs;
            end
            
            if curr_howling_time - last_gain_update_time < howling_detect_release_time/fs
                if disp_details_flag; disp(['No Howling Check ' num2str( curr_howling_time - last_gain_update_time ) ' seconds after gain change']); end
                flag_gain_control_coping_fail = 1;
            else
                %- num_of_frames_since_last_detection
                [num_of_frames_since_last_detection] = CalculateNumOfFramesFromLastDetection(curr_howling_time, last_gain_update_time, MSD, GC);
                
                %- control
                if num_of_frames_since_last_detection < MSD.Const.num_of_frames
                    num_of_frames_since_last_detection = max(num_of_frames_since_last_detection, MSD.FA_Detect.min_num_of_frames_since_last_detection);
                    
                    MSD_copy.magnitude_history_mat = MSD_copy.magnitude_history_mat(:, (end - num_of_frames_since_last_detection + 1):end);
                end
            end
        end
        
        if ~flag_gain_control_coping_fail
            [nBin, candidates, measurements] = MsdBasedHowlingDetect(MSD_copy, fs);
            
            % iii) Handle Candidates:
            %      * Handling the Gain-Control Mechanism according to the 'howling candidates'.
            % spectrum(m(1:n));
            % disp(min(meanSquared))
            if ~isempty(candidates)
                [noticeable_howling_candidates, noticeable_howling_measurements] = ...
                    ExamineHowlingCandidates(nBin, candidates, measurements,...
                                             n, fs,...
                                             MSD, GC,...
                                             disp_details_flag, is_plot_flag);
                % noticeable_howling_candidates = candidates;
                % noticeable_howling_measurements = measurements;
                
                if ~isempty(noticeable_howling_candidates)
                    if disp_details_flag
                        mean_energey_per_candidate = measurements(:, 2);
                        disp(['max(mean_energey_per_candidate) = ' num2str(max(mean_energey_per_candidate))]);
                    end
                    GC.LevelDetect.gain_dest = GC.LevelDetect.gain_dest/2;
                    
                    GC.LevelDetect.gain_max = GC.LevelDetect.gain_max - 0.5;
                    if GC.LevelDetect.gain_max < GC.LevelDetect.gain_max_lower_boundary
                        GC.LevelDetect.gain_max = GC.LevelDetect.gain_max_lower_boundary;
                    end
                    GC.Time.n_howling = n;
                    GC.Time.n_next_renew = GC.Time.n_howling + gain_control_renew_time;
                    
                    % HowlingFreqList:
                    HowlingFreqList = [HowlingFreqList;...
                        [ repmat( (n-1)/fs , length(noticeable_howling_candidates), 1),...
                        noticeable_howling_measurements]];
                end
            end
        end
    end
end

%% 3) Gain Smoothing:
if n >= GC.LevelDetect.level_detection_start_time
    %- reversed smoothed-decoupled peak-detector:
    [gain_smoothed, gain_r] = smoothedDecoupledPeakDetector(...
        -GC.LevelDetect.gain_dest, -GC.LevelDetect.gain_curr, -GC.LevelDetect.gain_r_old,...
        gain_alpha_a, gain_alpha_r);
    GC.LevelDetect.gain_curr = -gain_smoothed;
    GC.LevelDetect.gain_r_old = -gain_r;
end

end

