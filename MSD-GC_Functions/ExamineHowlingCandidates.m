function [noticeable_howling_candidates, noticeable_howling_measurements] = ...
    ExamineHowlingCandidates(nBin, candidates, measurements,...
    n, fs,...
    MSD, GC,...
    disp_details_flag, is_plot_flag)
% This function examines the 'howling candidates' with respect to:
% - Whether they cross a loudness threshold -
%     so they are not noise and about to be noticeable by the human-ear.
% - The status of the Gain-Control Mechanism.
% 
% Input:
% - nBin - all possible frequency bins (if peak-picking not active, algorithms search all bins).
% - candidates - indices of suspected frequency-howl bins.
% - measurements - of candidate frequency-howls.
% #
% - n - sample index
% - fs - sampling frequency
% #
% - MSD - struct
% - GC - struct
% #
% - disp_details_flag - display details if disp_details_flag > 0
% - is_plot_flag - plot figures for debugging if is_plot_flag > 1
% 
% Output:
% - noticeable_howling_candidates - which are not noise and about to be noticeable by the human-ear.
% - noticeable_howling_measurements - of "noticeable" candidate frequency-howls.
% 
% Author: Yehav Alkaher.
%%
if nargin < 8
    disp_details_flag = 0;
end
if nargin < 9
    is_plot_flag = 0;
end
%% Parsing
msd_frame_length = MSD.Const.msd_frame_length;
minimal_howl_energy_threshold = MSD.Const.minimal_howl_energy_threshold;
is_single_minimal_howl_energy_threshold = (length(minimal_howl_energy_threshold(:)) == 1);
is_per_freq_minimal_howl_energy_threshold = (length(minimal_howl_energy_threshold(:)) == msd_frame_length/2);
assert(is_single_minimal_howl_energy_threshold | is_per_freq_minimal_howl_energy_threshold)

howling_detect_release_time = GC.Time.Const.howling_detect_release_time;
%% Handling
if disp_details_flag
    disp('candidate howling freqs:')
    disp(num2str((candidates-1)*fs/msd_frame_length));
    disp(['time of suspecteed howling event: ' num2str(n/fs)])
end

% Remove unnoticeable howls:
mean_energey_per_candidate = measurements(:, 2);
if is_single_minimal_howl_energy_threshold
    noticeable_howling_candidates = candidates(mean_energey_per_candidate >= minimal_howl_energy_threshold);% dB
    noticeable_howling_measurements = measurements(mean_energey_per_candidate >= minimal_howl_energy_threshold, :);
else
    noticeable_howling_candidates = candidates(mean_energey_per_candidate >= minimal_howl_energy_threshold(candidates));% dB
    noticeable_howling_measurements = measurements(mean_energey_per_candidate >= minimal_howl_energy_threshold(candidates), :);
end
% noticeable_howling_candiates = candidates;
% noticeable_howling_measurements = measurements;

if isempty(noticeable_howling_candidates)
    noticeable_howling_measurements = [];
    if (is_plot_flag > 1)
        normalized_psd_zero2pi_db = CalculateNormalizedPSDdB(MSD.msd_buffer + db2mag(-600)*randn(size(MSD.msd_buffer)));
        zero2piRange = 1:length(normalized_psd_zero2pi_db);
        figure; plot((zero2piRange-1)*fs/msd_frame_length, normalized_psd_zero2pi_db)
        close
    end
    
    if disp_details_flag
        if is_single_minimal_howl_energy_threshold
            warning(['max(mean_energey_per_candidate) = ' num2str(max(mean_energey_per_candidate)) '< ' num2str(minimal_howl_energy_threshold) ' dB'])
        else
            warning(['max(mean_energey_per_candidate) = ' num2str(max(mean_energey_per_candidate)) 'dB'])
        end
    end
elseif n < GC.Time.n_howling + howling_detect_release_time
    if disp_details_flag; warning(['Howling still appears after' num2str( (n - GC.Time.n_howling)/fs ) 'seconds']); end
    noticeable_howling_candidates = [];
    noticeable_howling_measurements = [];
elseif all( (noticeable_howling_candidates-1)*fs/msd_frame_length < 15)% Hz
    [~, iValid, ~] = intersect(nBin, noticeable_howling_candidates);
    if (is_plot_flag > 1)
        figure;
        imagesc(MSD.magnitude_history_mat);
        ax = gca; ax.YDir = 'normal';
        hold on;
        plot(repelem(MSD.Const.num_of_frames, length(noticeable_howling_candidates)),...
             nBin(iValid), 'ro', 'MarkerSize', 10);
        hold off;
        close
    end
    if disp_details_flag; warning('All noticeable frequency-howls are below 15 Hz'); end
    noticeable_howling_candidates = [];
    noticeable_howling_measurements = [];
end
end
