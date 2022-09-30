%% Synopsis:
% * This script demonstrates the advantage of using two microphones
%   for the purpose of 'Adaptive Speech Reinforcement for In-Car Communications'.
% * Microphone Objectives:
%   - The first microphone is responsible for acquiring the speaker's speech (the driver)
%     and reproducing it via the loudspeaker.
%   - The second microphone is responsible for detecting the howling effect
%     (and even detecting positive feedback).
% * A note regarding the MSG (Maximal Stable Gain):
%   - The MSG is dependent in the RIR and not in the microphone placement.
%   - The 2nd microphone's objective is to better control the amplification gain.
%   - The 2nd microphone is supposed to increase the sensitivity for howling detection.
%   - Our objective is to show that the output speech signal reaches the desired gain over time.
%
% * The scenario in the script:
%   - 'Speaker' inside a 'Room' with a certain 'Environmental Noise'.
%   - A desired 'Amplification Gain' is chosen.
%
% * The room is tested for:
%   - Single Microphone.
%   - Two Microphones - where the 2nd microphone is somewhere in the room.
%
% * Desired Results:
%   - 'Hearing' the advantage in performance, for our solution and for the standard practice.
%   - Our objective is to show that the output speech signal reaches the desired gain over time.
%
% Author: Yehav Alkaher.

%%
clc; clear; close all;
addpaths_script
format shortg
%% 1. Choose Method:
method_flag = 'MSD_GainControl';
% method_flag = 'SimpleGain';
%% 2. Scenario Setting:
%% 2.1. Initializations:
Speech_dir_path = '';
Speech_filename = '564352__anzbot__ladies-and-gentlemen.mp3';

RIRs_dir_path = 'meshgrid_RIRs_ToTest\';
RIRs = load('meshgrid_RIRs_ToTest\generatedRIR_CarCabinComm_L_231.mat');

targetFolder = 'Results\Demo\';
if ~isdir(targetFolder)
    mkdir(targetFolder)
end

targetFolder = [targetFolder method_flag '\'];
if ~isdir(targetFolder)
    mkdir(targetFolder)
end

% Initialization Script:
Initialize_Speech_Reinforcement_System;

%% 2.2. Operation Configurations:
Conf.EnableNoise = 1;

Conf.isPlotFlag = 1;
Conf.isDispFlag = 0;%1;

Conf.Delays.systemTimeDelay = 0.005;%[Sec]
Conf.Delays.minPropagationTimeDelay = 0.003;%[Sec]

%% 2.3. Value Ranges:
single_mic_position_index = find((RIRs.mic2_X_grid(:) == 0.375) & (RIRs.mic2_Y_grid(:) == 2.500));

worst_case_pf_index_1     = find((RIRs.mic2_X_grid(:) == 0.125) & (RIRs.mic2_Y_grid(:) == 2.750));
worst_case_pf_index_2     = find((RIRs.mic2_X_grid(:) == 0.625) & (RIRs.mic2_Y_grid(:) == 0.500));
worst_case_pf_index_3     = find((RIRs.mic2_X_grid(:) == 0.625) & (RIRs.mic2_Y_grid(:) == 2.500));
worst_case_pf_index_4     = find((RIRs.mic2_X_grid(:) == 0.875) & (RIRs.mic2_Y_grid(:) == 1.000));
worst_case_pf_index_5     = find((RIRs.mic2_X_grid(:) == 1.125) & (RIRs.mic2_Y_grid(:) == 1.750));
worst_case_pf_index_6     = find((RIRs.mic2_X_grid(:) == 1.375) & (RIRs.mic2_Y_grid(:) == 1.500));
worst_case_pf_index_7     = find((RIRs.mic2_X_grid(:) == 1.625) & (RIRs.mic2_Y_grid(:) == 1.000));
worst_case_pf_index_8     = find((RIRs.mic2_X_grid(:) == 1.875) & (RIRs.mic2_Y_grid(:) == 0.750));


RIRsList = [single_mic_position_index,...
    worst_case_pf_index_1,...
    worst_case_pf_index_2,...
    worst_case_pf_index_3,...
    worst_case_pf_index_4,...
    worst_case_pf_index_5,...
    worst_case_pf_index_6,...
    worst_case_pf_index_7,...
    worst_case_pf_index_8];
kList = 6;%[1 3 6 7];
SnrList = 200;%[200, 40];% [dB]
%% Assesing Reference Gains for 'MSG Estimation':
% method_flag = 'SimpleGain'; RIRsList = single_mic_position_index; kList = [1 5 6]; SnrList = 200;
% Simulating Positive Feedback:
% method_flag = 'SimpleGain'; RIRsList = single_mic_position_index; kList = 8; SnrList = 200;

% RIRsList = [single_mic_position_index,...
%     worst_case_pf_index_1, worst_case_pf_index_2,...
%     worst_case_pf_index_3, worst_case_pf_index_4];
% kList = [1 3 4 5 6 7 8];
% SnrList = [200 40];

%%
all_measurements = [];
all_measurements_scenarios = {};

saveFiguresFlag = 1;

%% 2.4. Read samples:
[y,Fs] = audioread([Speech_dir_path Speech_filename]);
y = y(:,1);
if (fs ~= Fs)
    [U,D] = rat(fs/Fs);
    signal = resample(y,U,D);
else
    signal = y;
end
signal = [1*db2mag(-200)*randn(1,1.5*fs), signal(:).'];
signal = signal + 0*db2mag(-220)*randn(size(signal));
% Simulating Noise + Silence: 
% signal = [db2mag(-40)*randn(1,0.5*fs), db2mag(-200)*randn(1,9.5*fs)];
% signal = signal(1:7*fs);

%% 3. Test
% Read Room Impulse Response:
for room_idx = RIRsList
    
    Extract_meshgrid_RIRs_To_Test_script
    
    %% Initializations before SIMULATIONS:
    % FIR filters:
    g1 = h_a_u1(:).';% Loudspeaker to Speaker's Mic (Mic1)
    g2 = h_a_u2(:).';% Loudspeaker to Howling Detection Mic (Mic2)
    
    g1 = conv(g1,Conf.remove_DC_filter,'same');
    g2 = conv(g2,Conf.remove_DC_filter,'same');
    
    u1 = conv(signal,h_s_u1);
    u2 = conv(signal,h_s_u2);
    filtersLength = length(h_a_u1);
    
    clear signal_echo rir h_s_u1 h_s_u2 h_a_u1 h_a_u2
    
    for k_idx = 1:length(kList)
        k = kList(k_idx);
        
        %% SNR:
        for snr_idx = 1:length(SnrList)
            SNR = SnrList(snr_idx);
            noiseStd = db2mag(-SNR)*std(u1);
            
            disp(['k = ' num2str(k) ' , SNR = ' num2str(SNR)]);
            all_measurements_scenarios = [all_measurements_scenarios ; {['k = ' num2str(k) ' , SNR = ' num2str(SNR)]}];
            
            %% Handling the echo effect:
            tic;
            [m1, m2, y, measurements, HowlingFreqTable] = ModelBasedSystemFuncByName(method_flag, Conf, u1, u2, k, g1, g2, filtersLength, noiseStd, noise_rng_seed);
            
            end_runtime = toc;
            disp([num2str(floor(end_runtime/60)) 'minutes and ' num2str(mod(end_runtime,60), 4) ' seconds'])
            
            %% Save Results:
            if (saveFiguresFlag == 1)
                scenario_dir = ['K_' num2str(k) '_SNR_' num2str(SNR) '\'];
                if ~isdir([targetFolder scenario_dir])
                    mkdir([targetFolder scenario_dir])
                end
                
                if contains(method_flag,'MSD')
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_short_Gain_control_report_wrt_spectrograms.eps'],'epsc')
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_short_Gain_control_report_wrt_spectrograms.png'],'png')
                    savefig(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_short_Gain_control_report_wrt_spectrograms.fig'])
                    close;
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Gain_control_report_wrt_spectrograms.eps'],'epsc')
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Gain_control_report_wrt_spectrograms.png'],'png')
                    close;
                end
                saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Spectrogram_and_frequency_domain_comparison.eps'],'epsc')
                saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Spectrogram_and_frequency_domain_comparison.png'],'png')
                close;
                figure; spectrum(y); colormap autumn
                saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Spectrogram.eps'],'epsc')
                saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_measurements_Spectrogram.png'],'png')
                close;
                [Retro_HowlingFreqTable, ~] = RetrospectiveHowlingDetectionFunc(Conf, m2, u1, y, 1);
                if ~isempty(Retro_HowlingFreqTable)
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_retrospective_howling_detection.eps'],'epsc')
                    saveas(gcf,[targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)) '_retrospective_howling_detection.png'],'png')
                end
                close all;
                
                audiowrite([targetFolder scenario_dir Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_' num2str(RIRs.mic2_X_grid(room_idx)) '_' num2str(RIRs.mic2_Y_grid(room_idx)),...
                    '_loudspeaker_signal_PESQ_' num2str(measurements.PESQ),...
                    '_IsDist' num2str(measurements.OverallSymmISDist),...
                    '_AvgIsDist' num2str(measurements.ShortTermSymmISDistAvg),...
                    '.wav'], y, Fs);
                
                all_measurements = [all_measurements, measurements];
                
                % save([targetFolder Speech_filename(1:end-4) '_' method_flag '_' num2str(k) '_' num2str(SNR) '_measurements.mat'], 'measurements');
            end
        end
    end
end

save([targetFolder 'workspace.mat']);
%% For Posterior Analysis:
% ^Choose specific 'relevant mic2-locations' out of the collected data.
extra_analysis_flag = 0;
% load('Results\Demo\MSD_GainControl\workspace.mat')
% extra_analysis_flag = 1;
%% Performance-Summary:
% See 'effect_of_speech_reinforcement_on_measures_test.m'
if (saveFiguresFlag == 1)
    mic2_locations_mat = [...
        repelem('( ',length(RIRsList),1),...
        num2str(RIRs.mic2_X_grid(RIRsList).'),...
        repelem(' , ',length(RIRsList),1),...
        num2str(RIRs.mic2_Y_grid(RIRsList).'),...
        repelem(' )', length(RIRsList),1)];
    mic2_locations_cell = mat2cell(mic2_locations_mat, ones(length(RIRsList),1), size(mic2_locations_mat,2));
    
    if extra_analysis_flag
        relevant_mic2_location_indices = find(ismember(mic2_locations_mat, ...
            {'( 0.375 ,  2.5 )',...
            '( 0.625 ,  2.5 )',...
            '( 0.875 ,    1 )',...
            '( 1.125 , 1.75 )',...
            }));
        mic2_locations_mat = mic2_locations_mat(relevant_mic2_location_indices, :);
        mic2_locations_cell = mic2_locations_cell(relevant_mic2_location_indices);
    end
    %% Mic2-Location Comparison per Scenario:
    for k_idx = 1:length(kList)
        k = kList(k_idx);
        for snr_idx = 1:length(SnrList)
            SNR = SnrList(snr_idx);
            disp(['k = ' num2str(k) ' , SNR = ' num2str(SNR)]);
            
            relevant_scenario_indices = (0:(length(RIRsList) - 1))*length(kList)*length(SnrList) + (k_idx - 1)*length(SnrList) + snr_idx;
            
            relevant_measurements = all_measurements(relevant_scenario_indices);
            curr_scenario = unique(all_measurements_scenarios(relevant_scenario_indices));
            if (length(curr_scenario) ~= 1) || ~strcmp(curr_scenario, ['k = ' num2str(k) ' , SNR = ' num2str(SNR)])
                error('wrong relevant_room_indices')
            end
            
            if extra_analysis_flag
                relevant_measurements = relevant_measurements(relevant_mic2_location_indices);
            end
            
            measurements_mat = table2array(struct2table(relevant_measurements)).';
            
            normalized_measurements_mat = measurements_mat;
            normalized_measurements_mat(2,:) = normalized_measurements_mat(2,:)/4.5;
            normalized_measurements_mat(3:end-1,:) = exp(-normalized_measurements_mat(3:end-1,:));
            normalized_measurements_mat(end,:) = normalized_measurements_mat(end,:)/k;
            
            figure;
            set(gcf,'Position', [100 100 1000 800])% set(gcf,'Position', get(0, 'Screensize'))
            bar(measurements_mat);
            xticklabels(fields(measurements))
            legend(mic2_locations_cell)
            set(gca, 'XTickLabelRotation', 90)
            ylim([0, k+0.8])
            saveas(gcf,[targetFolder 'DemoSummaryGraph_Complete.png'],'png')
            saveas(gcf,[targetFolder 'DemoSummaryGraph_Complete.eps'],'epsc')
            
            figure;
            set(gcf,'Position', [100 100 1000 800])
            bar(normalized_measurements_mat);
            xticklabels(fields(measurements))
            legend(mic2_locations_cell)
            set(gca, 'XTickLabelRotation', 90)
            ylim([0 1.2])
            saveas(gcf,[targetFolder 'DemoSummaryGraph_Normalized.png'],'png')
            saveas(gcf,[targetFolder 'DemoSummaryGraph_Normalized.eps'],'epsc')
            
            %% Specific Feilds:
            figure;
            set(gcf,'Position', [50 50 1280 900])% set(gcf,'Position', get(0, 'Screensize'))
            subplot(3,1,1); bar(diag([relevant_measurements.PESQ]), 'stacked');
            xticklabels(mic2_locations_cell)
            xlabel('PESQ')
            ylim([0, 5]); grid on; grid minor;
            title('Quality - Maximize')
            subplot(3,1,2);
            bar([...
                [relevant_measurements.ShortTermISDistMedian];...
                [relevant_measurements.ShortTermISMeasMedian];...
                [relevant_measurements.ShortTermLLRMeasMedian];...
                [relevant_measurements.STFT_normalized_distorion];...
                [relevant_measurements.MSE_normalized_distorion]]);
            xticklabels({...
                'Short-Term\newlineIS-Dist Median',...
                'Short-Term\newlineIS-Meas Median',...
                'Short-Term\newlineLLR Median',...
                'STFT\newlineNormalized-Distorion',...
                'MSE\newlineNormalized-Distorion'})
            grid on; grid minor;
            title(['Distortion - Minimize'])
            legend(mic2_locations_cell)
            subplot(3,1,3);
            bar([...
                [relevant_measurements.ShortTermEffectiveGainMedian];...
                [relevant_measurements.STFT_matched_effective_gain];...
                [relevant_measurements.MSE_matched_effective_gain]]);
            xticklabels({...
                'Short-Term\newlineEffective-Gain Median',...
                'STFT-based\newlineMatched-Effective-Gain',...
                'MSE-based\newlineMatched-Effective-Gain'})
            grid on; grid minor;
            title(['Gain = ' num2str(k) ' - Maximize'])
            legend(mic2_locations_cell)
            
            suptitle(['Gain = ' num2str(k) ' , SNR = ' num2str(SNR) ' [dB]']);
            
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Selected.png'],'png')
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Selected.eps'],'epsc')
            
            %% Chosen Feilds:
            figure;
            set(gcf,'Position', [50 50 1280 900])% set(gcf,'Position', get(0, 'Screensize'))
            subplot(3,1,1); bar(diag([relevant_measurements.PESQ]), 'stacked');
            xticklabels(mic2_locations_cell)
            % xlabel('PESQ')
            ylim([0, 5]); grid on; grid minor;
            title(['Quality - Maximize\newline'...
                '*PESQ = [' num2str(round([relevant_measurements.PESQ], 3)) '] *'])
            subplot(3,1,2);
            bar(diag([relevant_measurements.ShortTermISDistMedian]), 'stacked');
            xticklabels(mic2_locations_cell)
            % xlabel('Short-Term IS-Dist Median')
            grid on; grid minor;
            title(['Distortion - Minimize\newline'...
                '*Short-Term IS-Dist Median = [' num2str(round([relevant_measurements.ShortTermISDistMedian], 3)) '] *'])
            subplot(3,1,3);
            bar(diag([relevant_measurements.MSE_matched_effective_gain]), 'stacked');
            xticklabels(mic2_locations_cell)
            % xlabel('MSE-based Matched-Effective-Gain')
            grid on; grid minor;
            title(['Gain = ' num2str(k) ' - Maximize\newline'...
                '*MSE-based Matched-Effective-Gain = [' num2str([relevant_measurements.MSE_matched_effective_gain]) '] *'])
            
            suptitle(['Gain = ' num2str(k) ' , SNR = ' num2str(SNR) ' [dB]']);
            
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Selected.png'],'png')
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Selected.eps'],'epsc')
            %% Final Chosen Feilds:
            figure;
            % set(gcf,'Position', [50 50 1280 900])% set(gcf,'Position', get(0, 'Screensize'))
            ax1 = subplot(2,1,1);
            bar(diag([relevant_measurements.ShortTermISDistMedian]), 'stacked');
            text(1:length([relevant_measurements.ShortTermISDistMedian]),[relevant_measurements.ShortTermISDistMedian],num2str([relevant_measurements.ShortTermISDistMedian].'),'vert','bottom','horiz','center'); 
            ylim([0, 1.15*max([relevant_measurements.ShortTermISDistMedian])])
            xticklabels(mic2_locations_cell)
            xlim([0.5, max(xticks)+0.5])
            % xlabel('Short-Term IS-Dist Median')
            grid on; grid minor;
            ylabel('Med-IS-Dist','interpreter','latex')
            ax2 = subplot(2,1,2);
            bar(diag([relevant_measurements.MSE_matched_effective_gain]), 'stacked');
            text(1:length([relevant_measurements.MSE_matched_effective_gain]),[relevant_measurements.MSE_matched_effective_gain],num2str([relevant_measurements.MSE_matched_effective_gain].'),'vert','bottom','horiz','center'); 
            ylim([0, 1.15*max([relevant_measurements.MSE_matched_effective_gain])])
            xticklabels(mic2_locations_cell)
            xlim([0.5, max(xticks)+0.5])
            % xlabel('MSE-based Matched-Effective-Gain')
            grid on; grid minor;
            % xlabel('MSE-based Matched-Effective-Gain $K_{{eff}}$','interpreter','latex')
            ylabel('$K_{{eff}}$','interpreter','latex')
            % suptitle(['Gain = ' num2str(k) ' , SNR = ' num2str(SNR) ' [dB]']);
            
            p1 = get(ax1, 'Position');
            p2 = get(ax2, 'Position');
            p1(2) = p2(2)+p2(4);
            set(ax1, 'pos', p1);
            
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Final_Latex.png'],'png')
            saveas(gcf,[targetFolder num2str(k) '_' num2str(SNR) '_DemoSummaryGraph_Final_Latex.eps'],'epsc')
        end
    end
    
    %% Amplification-Gain & Mic2-Location Comparison per SNR:
    for snr_idx = 1:length(SnrList)
        SNR = SnrList(snr_idx);
        disp(['SNR = ' num2str(SNR)]);

        k_idx = 1:length(kList);
        relevant_scenario_indices = (0:(length(RIRsList) - 1))*length(kList)*length(SnrList) + (k_idx.' - 1)*length(SnrList) + snr_idx;

        relevant_measurements = all_measurements(relevant_scenario_indices);
        curr_scenarios = unique(all_measurements_scenarios(relevant_scenario_indices));
        if length(curr_scenarios) ~= length(kList)
            error('wrong relevant_room_indices')
        end

        %% For Latex:
        figure;
        set(gcf,'Position', [100 50 600 870])% set(gcf,'Position', get(0, 'Screensize'))
        subplot(2,1,1);
        Z = arrayfun(@(x) x.ShortTermISDistMedian, relevant_measurements);
        b = bar3(Z);
        colorbar
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        xticklabels(mic2_locations_cell)
        xlabel('mic2-location')
        yticklabels(num2cell(kList))
        ylabel('Gain')
        zlabel('Med-IS-Dist','interpreter','latex')
        zlim([0, 1.15*max([all_measurements.ShortTermISDistMedian])])
        
        subplot(2,1,2);
        Z = arrayfun(@(x) x.MSE_matched_effective_gain, relevant_measurements);
        b = bar3(Z);
        colorbar
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        xticklabels(mic2_locations_cell)
        xlabel('mic2-location')
        yticklabels(num2cell(kList))
        ylabel('Gain')
        zlabel('$K_{{eff}}$','interpreter','latex')
        zlim([0, 1.15*max([all_measurements.MSE_matched_effective_gain])])
        
        suptitle(['Gain = [' num2str(kList) '] , SNR = ' num2str(SNR) ' [dB]']);

        saveas(gcf,[targetFolder 'SNR_' num2str(SNR) '_Combined_DemoSummaryGraph_Final_Latex.png'],'png')
        saveas(gcf,[targetFolder 'SNR_' num2str(SNR) '_Combined_DemoSummaryGraph_Final_Latex.eps'],'epsc')
        savefig(gcf,[targetFolder 'SNR_' num2str(SNR) '_Combined_DemoSummaryGraph_Final_Latex.fig'])
    end
end