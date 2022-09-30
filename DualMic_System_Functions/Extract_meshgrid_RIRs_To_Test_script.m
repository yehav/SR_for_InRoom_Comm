% This script extracts the room inpulse responses, created in
% 'RirGeneratorTest_CarCabinComm_Main1', for the Howling Reinforcement System scripts.

% Inputs:
% - RIRs_dir_path
% - RIRs_filenames
% - room_idx
% - fs

% Author: Yehav Alkaher.

%%
[row, col] = ind2sub(RIRs.mic2_grid_size,room_idx);
mic2_loc_subscripts = [row, col];
mic2_loc_x_y = [RIRs.mic2_X_grid(room_idx), RIRs.mic2_Y_grid(room_idx)];
mic2_loc_str = ['[' num2str(mic2_loc_x_y(1),'%.3f ') ' ' num2str(mic2_loc_x_y(2),'%.3f ') ' 1.000]'];
disp(['mic2_loc = ' mic2_loc_str])
Fs = RIRs.fs;
h_s_u1 = RIRs.h_s_u1;% Source to the System Receiver.
h_s_u2 = RIRs.source_to_mic2_RIR_cell{room_idx};% Source to Howling Detection Receiver.
h_s_u3 = RIRs.h_a_u3;% Source to Listener.

h_a_u1 = RIRs.h_a_u1;% Amplifier to the System Receiver.
h_a_u2 = RIRs.amp_to_mic2_RIR_cell{room_idx};% Amplifier to Howling Detection Receiver.
h_a_u3 = RIRs.h_a_u3;% Amplifier to Listener.

if (fs ~= Fs)
    [U,D] = rat(fs/Fs);
    h_s_u1 = resample(h_s_u1,U,D);
    h_s_u2 = resample(h_s_u2,U,D);
    h_s_u3 = resample(h_s_u3,U,D);
    
    h_a_u1 = resample(h_a_u1,U,D);
    h_a_u2 = resample(h_a_u2,U,D);
    h_a_u3 = resample(h_a_u3,U,D);
end

clear RIR
