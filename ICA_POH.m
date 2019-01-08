function [BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF)
% ICA_POH The Independent Component Analysis (ICA) Method from: Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated cardiac pulse measurements using video imaging and blind source separation. Optics express, 18(10), 10762-10774. DOI: 10.1364/OE.18.010762
%
%   Inputs:
%       VideoFile               = Video file path.
%       FS                      = Video framerate (fps).
%       StartTime               = Timepoint at which to start process (default = 0 seconds).
%       Duration                = Duration of the time window to process (default = 60 seconds).
%       ECGFile                 = File path to corresponding ECG data file (.mat) containing: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.
%       PPGFile                 = File path to corresponding PPG data file (.mat) containing: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       BVP                     = Processed Blood Volume Pulse (BVP).
%       PR                      = Estimated Pulse Rate (PR) from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate (HR) measured from the ECG timeseries R-waves for the window.
%       PR_PPG                  = Pulse Rate (PR) measured from the PPG timeseries systolic onsets for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio (SNR) calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013.
%
%   Requires - Signal Processing Toolbox
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, December 2018
% Copyright (c)
% Licensed under the MIT license.

addpath(genpath('tools'))

%% Parameters
LPF = 0.7; %low cutoff frequency (Hz) - 0.7 Hz in reference
HPF = 2.5; %high cutoff frequency (Hz) - 4.0 Hz in reference

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
T = zeros(FramesToRead,1);%initialize time vector
RGB=zeros(FramesToRead,3);%initialize color signal
FN=0;
while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in reference as using Viola Jones, 2004
    VidROI = VidFrame;
    
    %position for optional skin segmentation
    
    RGB(FN,:) = sum(sum(VidROI));%if different size regions are used for different frames, the signals should be normalized by the region size, but not necessary for whole frame processing or constant region size
end

%% Detrend & ICA:
NyquistF = 1/2*FS;
RGBNorm=zeros(size(RGB));
Lambda=100;
for c=1:3
    RGBDetrend= spdetrend(RGB(:,c),Lambda); %M. P. Tarvainen, TBME, 2002
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
load(ECGFile);
ECG.time = (1:length(ECG.data))/ECG.fs;
ECGMask = (ECG.time>=StartTime) & (ECG.time<=StartTime+Duration);
ECGPeakMask = ((ECG.peaks./ECG.fs)>=StartTime) & ((ECG.peaks./ECG.fs)<=StartTime+Duration);
HR_ECG = (1/mean(diff(ECG.peaks(ECGPeakMask)./ECG.fs)))*60;

if ~isempty(PPGFile)
    load(PPGFile);
    PPG.time = (1:length(PPG.data))/PPG.fs;
    PPGMask = (PPG.time>=StartTime) & (PPG.time<=StartTime+Duration);
    if isfield(PPG,'peaks')
        PPGPeakMask = ((PPG.peaks./PPG.fs)>=StartTime) & ((PPG.peaks./PPG.fs)<=StartTime+Duration);
        PR_PPG = (1/mean(diff(PPG.peaks(PPGPeakMask)./PPG.fs)))*60;
    else
        PR_PPG = NaN;
    end
else
    PR_PPG = NaN;
end

%% SNR
SNR = bvpsnr(BVP,FS,HR_ECG,PlotSNR);

%% Optionally Plot Timeseries
if(PlotTF)
    %Plot ECG, PPG, iPPG timeseries    
    figure
    
    if ~isempty(PPGFile)
        %Plot ECG
        Ax1=subplot(3,1,1);
        plot(ECG.time(ECGMask),ECG.data(ECGMask))
        hold on
        plot(ECG.peaks(ECGPeakMask)/ECG.fs,ECG.data(ECG.peaks(ECGPeakMask)),'*')
        ylabel('ECG (a.u.)')
        title('ICA Method - ECG, PPG, iPPG Timeseries')
        
        %Plot PPG
        Ax2=subplot(3,1,2);
        plot(PPG.time(PPGMask),PPG.data(PPGMask))
        hold on
        plot(PPG.peaks(PPGPeakMask)/PPG.fs,PPG.data(PPG.peaks(PPGPeakMask)),'*')
        ylabel('PPG (a.u.)')
        
        %Plot iPPG
        Ax3=subplot(3,1,3);
        plot(T,BVP)
        hold on
        ylabel('iPPG (a.u.)')

        xlabel('Time (s)')
        
        linkaxes([Ax1,Ax2,Ax3],'x')
    else
        %Plot ECG
        Ax1=subplot(2,1,1);
        plot(ECG.time(ECGMask),ECG.data(ECGMask))
        hold on
        plot(ECG.peaks(ECGPeakMask)/ECG.fs,ECG.data(ECG.peaks(ECGPeakMask)),'*')
        ylabel('ECG (a.u.)')
        title('ICA Method - ECG, iPPG Timeseries')
        
        %Plot iPPG
        Ax2=subplot(2,1,2);
        plot(T,BVP)      
        hold on
        ylabel('iPPG (a.u.)')

        xlabel('Time (s)')
        
        linkaxes([Ax1,Ax2],'x')
    end
    
    xlim([StartTime StartTime+Duration])

end%endif plot

end%end function
