function [BVP, PR, HR_ECG, PR_PPG, SNR] = CHROM_DEHAAN(VideoFile, FS, StartTime, Duration, ECGData, PPGData, PlotTF)
% CHROM_DEHAAN The CHROM (DeHaan et al. 2014) Method.
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
%       BVP                     = Processed Blood Volume Pulse using the chrominance method (DeHaan et al. 2014).
%       PR                      = Estimated Pulse Rate from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate measured from the ECG timeseries for the window.
%       PR_PPG                  = Pulse Rate measured from the PPG timeseries for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013
%
%   Requires - Signal Processing Toolbox
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, December 2018

%% Parameters
SkinSegmentTF=false;

LPF = 0.70; %low cutoff frequency (Hz)
HPF = 4.0; %high cutoff frequency (Hz)
% NOTE DIFFERENT FROM AS IN THE PAPER FOR CONSISTENCY:
%LPF = 0.667; %low cutoff frequency (Hz) - specified as 40 bpm
%HPF = 4.0; %high cutoff frequency (Hz) - specified as 240 bpm

WinSec=1.6;%(was a 32 frame window with 20 fps camera)

%% Add Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)
    %addpath([cd '\optional\']);%GNU GPL rgb2ycbcr.m function
    addpath([cd '\optional\rgb2ycbcr.m']);%GNU GPL rgb2ycbcr.m function
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

FramesToRead=floor(Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

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
    
    if(SkinSegmentTF)%skin segmentation - not specified in reference
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
    else
        RGB(FN,:) = sum(sum(VidROI,2))./(size(VidROI,1)*size(VidROI,2));
    end
end%endwhile video

if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
    rmpath([cd '\optional\']);
end

%% CHROM:
NyquistF = 1/2*FS;
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter - originally specified as an FIR filter

%Window parameters - overlap, add with 50% overlap
WinL = ceil(WinSec*FS);
if(mod(WinL,2))%force even window size for overlap, add of hanning windowed signals
    WinL=WinL+1;
end
NWin = floor((FN-WinL/2)/(WinL/2));
S = zeros(NWin,1);
WinS=1;%Window Start Index
WinM=WinS+WinL/2;%Window Middle Index
WinE=WinS+WinL-1;%Window End Index

for i = 1:NWin
    TWin = T(WinS:WinE,:);
    
    RGBBase = mean(RGB(WinS:WinE,:));
    RGBNorm = bsxfun(@times,RGB(WinS:WinE,:),1./RGBBase)-1;
    
    % CHROM
    Xs = squeeze(3*RGBNorm(:,1)-2*RGBNorm(:,2));%3Rn-2Gn
    Ys = squeeze(1.5*RGBNorm(:,1)+RGBNorm(:,2)-1.5*RGBNorm(:,3));%1.5Rn+Gn-1.5Bn
    
    Xf=filtfilt(B,A,double(Xs));
    Yf=filtfilt(B,A,double(Ys));
    
    Alpha=std(Xf)./std(Yf);
    
    SWin = Xf - Alpha.*Yf;
    
    SWin = hann(WinL).*SWin;
    %overlap, add Hanning windowed signals
    if(i==1)
        S=SWin;
        TX=TWin;
    else
        S(WinS:WinM-1)=S(WinS:WinM-1)+SWin(1:WinL/2);%1st half overlap
        S(WinM:WinE)=SWin(WinL/2+1:end);%2nd half
        TX(WinM:WinE)=TWin(WinL/2+1:end);
    end
    
    WinS=WinM;
    WinM=WinS+WinL/2;
    WinE=WinS+WinL-1;
end

BVP=S;
T=T(1:length(BVP));

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

%% Remove Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
    rmpath([cd '\optional\']);
end

end%end function
