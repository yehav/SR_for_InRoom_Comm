# Speech Reinforcement for In-Room Communications
## by Yehav Alkaher
## Preface
This code implements a speech reinforcement system & algorithm.
The files include two main demonstration scripts that correspond to papers:
* Temporal Howling Detector for Speech Reinforcement Systems
	* Authors:
	Yehav Alkaher, Prof. Israel Cohen
		* Andrew and Erna Viterbi Faculty of Electrical and Computer Engineering, Technion – Israel Institute of Technology, Haifa, Israel
* Dual-Microphone Speech Reinforcement System with Howling-Control for In-Car Speech Communication
	* Authors:
	Yehav Alkaher, Prof. Israel Cohen
		* Andrew and Erna Viterbi Faculty of Electrical and Computer Engineering, Technion – Israel Institute of Technology, Haifa, Israel
	* Citation:
Alkaher Y and Cohen I (2022) Dual-Microphone Speech Reinforcement System With Howling-Control for In-Car Speech Communication. _Front. Sig. Proc._ 2:819113. doi: 10.3389/frsip.2022.819113


## Code
MATLAB Version: R2017b
### Main scripts and functions:
* Temporal Howling Detector Demo:
	* 'effect_of_two_pole_TF_main_script.m'
	* 'MSD-GC_Functions'
* Dual-Microphone Speech-Reinforcement System & Algorithm:
	* 'dual_mic_MSD_GC_demonstration_script.m'
	* 'DualMic_System_Functions'
* Utils

#### Inner-Function References:
	* Room Impulse Responses (RIRs) extracted via the RIR 	Generator:
		Habets, E. A. (2006). Room impulse response generator. Technische Universiteit Eindhoven, Tech. Rep 2, 1
	* ISO226:2003:
		Christopher Hummersone (2021). ISO 226:2003 Normal equal-loudness-level contours (https://github.com/IoSR-Surrey/MatlabToolbox), GitHub. Retrieved July 11, 2021.
	* STOI:
		Copyright 2009: Delft University of Technology, Signal & Information Processing Lab. The software is free for non-commercial use. This program comes WITHOUT ANY WARRANTY.
	* PESQ:
		Jacob Donley 2017.

#### Sample Audio was downloaded from:
	Font F, Roma G, Serra X. Freesound technical demo. In: MM'13. Proceedings of the 21st ACM international conference on Multimedia; 2013 Oct 21-25; Barcelona, Spain. New York: ACM; 2013. p. 411-2. DOI: 10.1145/2502081.2502245

