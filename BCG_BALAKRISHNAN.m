function [BCG, PR, HR_ECG, PR_PPG, SNR] = BCG_BALAKRISHNAN(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% BCG_BALAKRISHNAN The BCG (Balakrishnan et al., 2012) Method applied to the AFRL dataset.
%
%   Inputs:
%       VideoFile               = Video filename.
%       StartTime               = Timepoint at which to start process (default = 15);
%       Duration                = Duration of the time window to process (default = 30 seconds).
%       BioSemiData             = Corresponding BioSemiData file.
%       ECGMark                 = Corresponding BioSemiData_ECGMark data file.
%       PPGMark                 = Corresponding BioSemiData_ECGMark data file.
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       BCG                     = Processed Balistocardiogram using (Balakrishnan et al., 2012) Method.
%       PR                      = Estimated Pulse Rate from processed BCG timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate measured from the ECG timeseries for the window.
%       PR_PPG                  = Pulse Rate measured from the PPG timeseries for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency using A method adapted from the method by G. de Haan, TBME, 2013
%
%   Requires - Computer Vision Toolbox
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, June 2018
%% Parameters
FS = 120; %true frame rate

LPF = 0.70; %low cutoff frequency (Hz)
HPF = 4.0; %high cutoff frequency (Hz)
% NOTE DIFFERENT FROM AS IN THE PAPER FOR CONSISTENCY:
%LPF = 0.75; %low cutoff frequency (Hz)
%HPF = 5; %high cutoff frequency (Hz)

%% Add Backup Functions
if(~license('test', 'Statistics_Toolbox'))
    addpath([cd '\optional\']);%GNU GPL quantile.m function
end

%% Plot Control
if(PlotTF)
    PlotPRPSD = true;
    PlotSNR = true;
else
    PlotPRPSD = false;
    PlotSNR = false;
end

%% Load Video:
VidObj = VideoReader(VideoFile);
VidObj.CurrentTime = StartTime;

FramesToRead=ceil(Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Initialize Tracking
BBox = [1,1,VidObj.width,VidObj.height];

MinQ = 0.01;
VidFrame = readFrame(VidObj);
Points0 = detectMinEigenFeatures(rgb2gray(VidFrame),'ROI',BBox,'MinQuality',MinQ);

if length(Points0) < 5
    BCG=NaN; PR=NaN; HR_ECG=NaN; PR_PPG=NaN; SNR=NaN;
    return
end

tracker = vision.PointTracker;
initialize(tracker,Points0.Location,VidFrame);

%% Read Video and Spatially Average:
T = zeros(FramesToRead-1,1);%Initialize Time Vector
D = zeros(FramesToRead-1,size(Points0.Location,1));%initialize D
FN=0;
while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in reference as using Viola Jones, 2004 to isolate face and region including the eyes
    
    [Points, Validity] = step(tracker,VidFrame);
    D(FN,:) = Points(:,2);%y, vertical component
end
%originally specified in reference with cubic-spine interpolation to 250Hz (from 30 Hz).

%% Signal Processing:
NyquistF = 1/2*FS;
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%5th order butterworth filter in reference
D_Filt = filtfilt(B,A,double(D));

D_Filt2 = bsxfun(@minus, D_Filt, D_Filt(1,:));
DL2 = sqrt(sum(D_Filt2.^2,2));
DMask = DL2<=quantile(DL2,.75);

D_Filt3=D_Filt2(DMask,:);
[Coeff, ~, Latent] = pca(D_Filt3);

Score = bsxfun(@minus, D_Filt, mean(D_Filt))/Coeff';

%% Component Selection
F=find(T>0);
[Pxx,F2] = plomb(Score(F,1:5),T(F));%only evaluate first 5 components

FMask = (F2 >= LPF)&(F2 <= HPF);
FRange = F2(FMask);
R = zeros(1,5);
for i = 1:5
    [MaxP,IDXP] = max(Pxx(:,i));
    R(i) = (MaxP+Pxx(IDXP*2,i))/sum(Pxx(:,i));
end
[~,NumPC] = max(R);

%% PR
PR = FRange(argmax(Pxx(FMask,NumPC),1))*60;

BCG_F = Score(:,NumPC);

BCG_N = BCG_F-mean(BCG_F); 
BCG=BCG_N;

%% Ground Truth HR:
ECG=load(ECGMark,'KeepBeatData');
%KeepBeatData Format (ECG)- Beat Index, Inter-Beat Interval (ms), Zero-Referenced Sample Index of R-wave Maximum, Zero-Referenced Time of R-Wave (s), Amplitude of R-Wave in Filtered Signal, Zero-Referenced Sample Index of R-wave Maximum, Zero-Referenced Time of S-Wave (s), Amplitude of S-Wave in Filtered Signal, User Added Data (0-automatic, 1-user added/ corrected)
ECGMask = (ECG.KeepBeatData(:,4)>=StartTime)&(ECG.KeepBeatData(:,4)<=StartTime+Duration);
HR_ECG = 1/mean(diff(ECG.KeepBeatData(ECGMask,4)))*60;

PPG=load(PPGMark,'KeepBeatData');
%KeepBeatData Format (PPG)- Pulse Index, Inter-Pulse Interval (ms), Zero-Referenced Sample Index of Systolic Onset, Zero-Referenced Time of Systolic Onset (s), Amplitude of Pulse/ Systolic Onset in Filtered Signal, Empty, Empty, Empty, User Added Data (0-automatic, 1-user added/ corrected)
PPGMask = (PPG.KeepBeatData(:,4)>=StartTime)&(PPG.KeepBeatData(:,4)<=StartTime+Duration);
PR_PPG = 1/mean(diff(PPG.KeepBeatData(PPGMask,4)))*60;

%% SNR
SNR = bvpsnr(BCG,FS,HR_ECG,PlotSNR);

%% Optionally Plot Timeseries
if(PlotTF)
    %Plot ECG, PPG, iPPG timeseries
    load(BioSemiData,'FiltNECG', 'FiltNPPG', 'ResampleTimeSeries')%load filtered and resampled ECG, PPG timeseries data
    
    figure
    
    Ax1=subplot(3,1,1);
    plot(ResampleTimeSeries,FiltNECG)
    hold on
    plot(ECG.KeepBeatData(:,4),ECG.KeepBeatData(:,5),'*')
    ylabel('ECG (a.u.)')
    title('ECG, PPG, iPPG Timeseries')
    
    Ax2=subplot(3,1,2);
    plot(ResampleTimeSeries,FiltNPPG)
    hold on
    plot(PPG.KeepBeatData(:,4),-1*PPG.KeepBeatData(:,5),'*')
    ylabel('PPG (a.u.)')
    
    Ax3=subplot(3,1,3);
    plot(T,BCG)
    hold on
    ylabel('BCG (a.u.)')
    
    xlabel('Time (s)')
    
    linkaxes([Ax1,Ax2,Ax3],'x')
    xlim([StartTime StartTime+Duration])
end%endif plot

%% Remove Backup Functions
if(~license('test', 'Statistics_Toolbox'))
    rmpath([cd '\optional\']);%GNU GPL quantile.m function
end

end%end function
