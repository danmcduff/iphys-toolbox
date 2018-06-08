function [BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% ICA_POH The ICA (Poh, McDuff, Picard 2010) Method applied to the AFRL dataset.
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
%       BVP                     = Processed Blood Volume Pulse using JADE ICA.
%       PR                      = Estimated Pulse Rate from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate measured from the ECG timeseries for the window.
%       PR_PPG                  = Pulse Rate measured from the PPG timeseries for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013
%
%   Requires - Signal Processing Toolbox
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, April 2018
%% Parameters
FS = 120; %true frame rate

LPF = 0.70; %low cutoff frequency (Hz)
HPF = 4.0; %high cutoff frequency (Hz)

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

%% Read Video and Spatially Average:
T = zeros(FramesToRead,1);%initialize Time Vector
RGB=zeros(FramesToRead,3);%initialize color signal
FN=0;
while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in reference as using Viola Jones, 2004
    VidROI = VidFrame;
    
    %position for optional skin segmentation
    
    RGB(FN,:) = sum(sum(VidROI));
end

%% Detrend & ICA:
NyquistF = 1/2*FS;
RGBNorm=zeros(size(RGB));
Lambda=1000;
for c=1:3
    RGBDetrend= spdetrend(RGB(:,c),Lambda);%M. P. Tarvainen, TBME, 2002
    RGBNorm(:,c) = (RGBDetrend - mean(RGBDetrend))/std(RGBDetrend); %normalize to zero mean and unit variance
end
[W,S] = ica(RGBNorm',3); %JADE ICA - J. F. Cardoso 1997, G. D. Clifford, MIT, 2004

%% Select BVP Source:
% Component with maximum normalized (by total power) power
MaxPx=zeros(1,3);
for c=1:3
    FF = fft(S(c,:));
    F=(1:length(FF))/length(FF)*FS*60;
    FF(1)=[];
    N=length(FF);
    Px = abs(FF(1:floor(N/2))).^2;
    Fx = (1:N/2)/(N/2)*NyquistF;
    Px=Px/sum(Px);
    MaxPx(c)=max(Px);
end

[M,MaxComp]=max(MaxPx(:));
BVP_I = S(MaxComp,:);

%% Filter, Normalize
%originally specified in reference with 5-point moving average, bandpass
%filter, and cubic-spine interpolation to 256Hz
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter
BVP_F = filtfilt(B,A,double(BVP_I));

BVP=BVP_F;

% Estimate Pulse Rate from periodogram
PR = prpsd(BVP,FS,40,240,PlotPRPSD);

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
SNR = bvpsnr(BVP,FS,HR_ECG,PlotSNR);

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
    plot(T,BVP)
    hold on
    ylabel('iPPG (a.u.)')
    
    xlabel('Time (s)')
    
    linkaxes([Ax1,Ax2,Ax3],'x')
    xlim([StartTime StartTime+Duration])
end%endif plot

end%end function