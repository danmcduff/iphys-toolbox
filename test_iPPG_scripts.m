clc, clear; close all
addpath(genpath([cd '\tools\']))%additional required functions

%% Set variables:
DataDirectory           = [cd '\test_data\'];%default determined off this script's directory
VideoFile               = [DataDirectory 'video_example.mp4'];% Video file path.
FS                      = 30;% Video framerate (fps).
StartTime               = 0;% Timepoint at which to start process (default = 0 seconds).
Duration                = 60;% Duration of the time window to process (default = 60 seconds).
ECGFile                 = [DataDirectory 'ECGData.mat'];% File path to corresponding ECG data file (.mat) containing: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.
PPGFile                 = [DataDirectory 'PPGData.mat'];% File path to corresponding PPG data file (.mat) containing: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.
PlotTF                  = true;% Boolean to turn plotting results on or off.

%% Check Files
if(~exist(VideoFile,'file'))
    error('VideoFile: ''%s'' does not exist.\nAn example video file can be downloaded from: ''https://drive.google.com/open?id=1oD4VbBD9ColSlbiIMEgxbvQ7LnXHPy1_''',VideoFile)
end
if(~exist(ECGFile,'file'))
    error('ECGFile: ''%s'' does not exist.\nAn example ECG file is contained in the respository ''test_data'' directory',ECGFile)
end
if(~exist(PPGFile,'file'))
    error('PPGGFile: ''%s'' does not exist.\nAn example PPG file is contained in the respository ''test_data'' directory',PPGFile)
end

%% Green - Verkruysse, Svaasand, Nelson (2008)
%Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote
%plethysmographic imaging using ambient light. Optics express, 16(26),
%21434-21445. DOI: 10.1364/OE.16.021434
[BVP, PR, HR_ECG, PR_PPG, SNR] = GREEN_VERKRUYSSE(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF);
fprintf('GREEN_VERKRUYSSE:\n')
PR, HR_ECG, PR_PPG, SNR
%% ICA - Poh, McDuff, Picard (2010)
% Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated
% cardiac pulse measurements using video imaging and blind source
% separation. Optics express, 18(10), 10762-10774. DOI: 10.1364/OE.18.010762
[BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF);
fprintf('ICA_POH:\n')
PR, HR_ECG, PR_PPG, SNR
%% CHROM - De Haan & Jeanne (2013)
% De Haan, G., & Jeanne, V. (2013). Robust pulse rate from
% chrominance-based rPPG. IEEE Transactions on Biomedical Engineering,
% 60(10), 2878-2886. DOI: 10.1109/TBME.2013.2266196
[BVP, PR, HR_ECG, PR_PPG, SNR] = CHROM_DEHAAN(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF);
fprintf('CHROM_DEHAAN:\n')
PR, HR_ECG, PR_PPG, SNR
%% POS - Wang, den Brinker, Stuijk & de Haan (2017)
% Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017).
% Algorithmic principles of remote PPG. IEEE Transactions on Biomedical
% Engineering, 64(7), 1479-1491. DOI: 10.1109/TBME.2016.2609282
[BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF);
fprintf('POS_WANG:\n')
PR, HR_ECG, PR_PPG, SNR
%% BCG - Balakrishnan, Durand & Guttag (2013)
% Balakrishnan, G., Durand, F., & Guttag, J. (2013, June). Detecting pulse
% from head motions in video. In Computer Vision and Pattern Recognition
% (CVPR), 2013 IEEE Conference on (pp. 3430-3437). IEEE. DOI: 10.1109/CVPR.2013.440
[BCG, PR, HR_ECG, PR_PPG, SNR] = BCG_BALAKRISHNAN(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF);
fprintf('BCG_BALAKRISHNAN:\n')
PR, HR_ECG, PR_PPG, SNR