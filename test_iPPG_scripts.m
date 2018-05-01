clear; close all; clc;

%% Set variables:
VideoFile               = 'E:\EX\test\P05T01VideoB2_MSMPEG4V3_Q95.avi';
StartTime               = 15;
Duration                = 30;
BioSemiData             = 'E:\EX\test\P05T01_BioSEMIData.mat';
ECGMark                 = 'E:\EX\test\P05T01_BioSEMIData_ECGMark.mat';
PPGMark                 = 'E:\EX\test\P05T01_BioSEMIData_PPGMark.mat';
PlotTF                  = true;

%% ICA - Poh, McDuff, Picard (2010)
% Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated
% cardiac pulse measurements using video imaging and blind source
% separation. Optics express, 18(10), 10762-10774.
[BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF);

%% CHROM - De Haan & Jeanne (2013)
% De Haan, G., & Jeanne, V. (2013). Robust pulse rate from
% chrominance-based rPPG. IEEE Transactions on Biomedical Engineering,
% 60(10), 2878-2886.
[BVP, PR, HR_ECG, PR_PPG, SNR] = CHROM_DEHAAN(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF);

%% POS - Wang, den Brinker, Stuijk & de Haan (2017)
% Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017).
% Algorithmic principles of remote PPG. IEEE Transactions on Biomedical
% Engineering, 64(7), 1479-1491.
[BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF);

%% BCG - Balakrishnan, Durand & Guttag (2013)
% Balakrishnan, G., Durand, F., & Guttag, J. (2013, June). Detecting pulse
% from head motions in video. In Computer Vision and Pattern Recognition
% (CVPR), 2013 IEEE Conference on (pp. 3430-3437). IEEE.
[BVP, PR, HR_ECG, PR_PPG, SNR] = BCG_BALAKRISHNAN(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF);
