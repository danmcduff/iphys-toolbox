clear; close all; clc;

VideoFile               = 'E:\EX\test\P05T01VideoB2_MSMPEG4V3_Q95.avi';
StartTime               = 15;
Duration                = 30;
BioSemiData             = 'E:\EX\test\P05T01_BioSEMIData.mat';
ECGMark                 = 'E:\EX\test\P05T01_BioSEMIData_ECGMark.mat';
PPGMark                 = 'E:\EX\test\P05T01_BioSEMIData_PPGMark.mat';
PlotTF                  = true;

[BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF);