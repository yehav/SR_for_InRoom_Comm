function [num_of_frames_since_last_detection] = CalculateNumOfFramesFromLastDetection(curr_howling_time, last_howling_time, MSD, GC)
% This is used in the process of Howling Detection False-Alarm Control.
% The function returns the num_of_frame since the last Howling-Detection,
% in order to consider only the relevant measurements,
% after the applied Gain is back with a large value.
% 
% Author: Yehav Alkaher.

fs = MSD.Const.fs;

num_of_frames_since_last_detection = inf;
if ~isempty(last_howling_time)
    time_diff_from_last_detection = curr_howling_time - last_howling_time;

    if time_diff_from_last_detection > GC.Time.Const.howling_detect_release_time/fs
        % disp(time_diff_from_last_detection);
        if time_diff_from_last_detection > GC.Time.Const.gain_control_renew_time/fs
            time_diff_from_last_gain_update = time_diff_from_last_detection - GC.Time.Const.gain_control_renew_time/fs;
        elseif time_diff_from_last_detection > GC.Time.Const.gain_control_release_time/fs
            time_diff_from_last_gain_update = time_diff_from_last_detection - GC.Time.Const.gain_control_release_time/fs;
        else
            time_diff_from_last_gain_update = time_diff_from_last_detection - GC.Time.Const.howling_detect_release_time/fs;
        end
        
        % 1 + (howling_history_delay*fs - MSD.Const.msd_frame_length) / MSD.Const.msd_frame_shift
        num_of_frames_since_last_detection = 1 + (time_diff_from_last_gain_update*fs - MSD.Const.msd_frame_length) / MSD.Const.msd_frame_shift;
    end
end

if num_of_frames_since_last_detection >= MSD.Const.num_of_frames
    num_of_frames_since_last_detection = MSD.Const.num_of_frames;
else
    num_of_frames_since_last_detection = floor(num_of_frames_since_last_detection);
end

end

