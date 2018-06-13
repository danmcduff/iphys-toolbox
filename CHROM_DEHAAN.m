function [BVP, PR, HR_ECG, PR_PPG, SNR] = CHROM_DEHAAN(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% CHROM_DEHAAN The CHROM (DeHaan et al. 2014) Method applied to the AFRL dataset.
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
%       BVP                     = Processed Blood Volume Pulse using the chrominance method (DeHaan et al. 2014).
%       PR                      = Estimated Pulse Rate from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate measured from the ECG timeseries for the window.
%       PR_PPG                  = Pulse Rate measured from the PPG timeseries for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013
%
%   Requires - Signal Processing Toolbox
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, June 2018

%% Parameters
FS = 120; %true frame rate

SkinSegmentTF=false;

LPF = 0.7; %low cutoff frequency (Hz) - specified as 40 bpm (0.667 Hz)
HPF = 2.5; %high cutoff frequency (Hz) - specified as 240 bpm (4.0 Hz)

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
    plot(TX,BVP)
    hold on
    ylabel('iPPG (a.u.)')
    
    xlabel('Time (s)')
    
    linkaxes([Ax1,Ax2,Ax3],'x')
    xlim([TX(1) TX(end)])%to account for differences caused by windows
end%endif plot

%% Remove Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
    rmpath([cd '\optional\']);
end

end%end function
