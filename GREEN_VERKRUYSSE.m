function [BVP, PR, HR_ECG, PR_PPG, SNR] = GREEN_VERKRUYSSE(VideoFile, FS, StartTime, Duration, ECGData, PPGData, PlotTF)
% GREEN_VERKRUYSSE The method (Verkruysse, Svaasand, Nelson 2008) Method.
% 
%   Inputs:
%       VideoFile               = Video filename.
%       FS                      = Video framerate.
%       StartTime               = Timepoint at which to start process (default = 15);
%       Duration                = Duration of the time window to process (default = 30 seconds).
%       ECGData                 = Corresponding ECGData data file.
%       PPGData                 = Corresponding PPGData data file.
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
% Daniel McDuff, Ethan Blackford, Justin Estepp, December 2018
% Copyright (c)
% Licensed under the MIT license.

addpath(genpath('tools'))

%% Parameters
LPF = 0.7; 
HPF = 2.5; 

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
    
    %position for optional face detection/tracking - originally specified in reference as a manual segmentation.
    VidROI = VidFrame;
    
    %position for optional skin segmentation
    
    RGB(FN,:) = sum(sum(VidROI));
end

%% Select BVP Source:
% Component with maximum normalized (by total power) power
BVP_I = RGB(:,2);

%% Filter, Normalize
%originally specified in reference with 5-point moving average, bandpass
%filter, and cubic-spine interpolation to 256Hz
NyquistF = 1/2*FS;
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter
BVP_F = filtfilt(B,A,double(BVP_I));

BVP=BVP_F;

% Estimate Pulse Rate from periodogram
PR = prpsd(BVP,FS,40,240,PlotPRPSD);

%% Ground Truth HR:
load(ECGData);
ECG.time = [1:length(ECG.data)]/ECG.fs;
ECGMask = (ECG.time>=StartTime)&(ECG.time<=StartTime+Duration);
ECGPeakMask = ((ECG.peaks./ECG.fs)>=StartTime)&((ECG.peaks./ECG.fs)<=StartTime+Duration);
HR_ECG = 1/mean(diff(ECG.peaks(ECGPeakMask)./ECG.fs))*60;

if ~isempty(PPGData)
    load(PPGData);
    PPG.time = [1:length(PPG.data)]/PPG.fs;
    PPGMask = (PPG.time>=StartTime)&(PPG.time<=StartTime+Duration);
    if isfield(PPG,'peaks')
        PPGPeakMask = ((PPG.peaks./PPG.fs)>=StartTime)&((PPG.peaks./PPG.fs)<=StartTime+Duration);
        PR_PPG = 1/mean(diff(PPG.data(PPGPeakMask)./PPG.fs))*60;
    else
        PR_PPG = nan;
    end
else
    PR_PPG = nan;
end

%% SNR
SNR = bvpsnr(BVP,FS,HR_ECG,PlotSNR);

%% Optionally Plot Timeseries
if(PlotTF)
    %Plot ECG, PPG, iPPG timeseries    
    figure
    
    if ~isempty(PPGData)
        %Plot ECG
        Ax1=subplot(3,1,1);
        plot(ECG.time(ECGMask),ECG.data(ECGMask))
        hold on
        plot(ECG.peaks(ECGPeakMask)/ECG.fs,ECG.data(ECG.peaks(ECGPeakMask)),'*')
        ylabel('ECG (a.u.)')
        title('ECG, PPG, iPPG Timeseries')
        %Plot PPG
        Ax2=subplot(3,1,2);
        plot(PPG.time(PPGMask),PPG.data(PPGMask))
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
        title('ECG, PPG, iPPG Timeseries')
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
