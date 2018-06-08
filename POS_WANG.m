function [BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% POS_WANG The POS (Wang et al. 2016) Method applied to the AFRL dataset.
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
%       BVP                     = Processed Blood Volume Pulse using Plane Orthogonal to Skin (POS) method (Wang et al. 2016).
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

LPF = 0.6667; %low cutoff frequency (Hz) - specified as 40 bpm (only used for pow
HPF = 4.0; %high cutoff frequency (Hz)

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
T = zeros(FramesToRead,1);%initialize Time Vector
RGB=zeros(FramesToRead,3);%initialize color signal
FN=0;
while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in
    %reference as a CSK detector from Henriques et al., 2012
    VidROI = VidFrame;
    
    if(SkinSegmentTF)%skin segmentation - originally specified in reference as an OC-SVM from Wang
    %et al. 2015
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

%% POS:
% Wenjin's transform
RGBBase = mean(RGB);
RGBNorm = bsxfun(@times,RGB,1./RGBBase)-1;
FF = fft(RGBNorm);
F = (0:size(RGBNorm,1)-1)*FS/size(RGBNorm,1);
S = FF*[-1/sqrt(6);2/sqrt(6);-1/sqrt(6)];
W = (S.*conj(S))./sum(FF.*conj(FF),2);
FMask = (F >= LPF)&(F <= HPF);%40-240 bpm
% FMask(length(FMask)/2+1:end)=FMask(length(FMask)/2:-1:1);
FMask = FMask + fliplr(FMask);
W=W.*FMask';%rectangular filter in frequency domain - not specified in original paper
FF = FF.*repmat(W,[1,3]);
RGBNorm=real(ifft(FF));
RGBNorm = bsxfun(@times,RGBNorm+1,RGBBase);

N = size(RGBNorm,1);
WinL = ceil(WinSec*FS);
S = zeros(N,1);
for i = 1:N-1
    if(i+WinL-1>length(S))%end of signal
        XWin = RGBNorm(i:length(S),:);
    else
        XWin = RGBNorm(i:i+WinL-1,:);
    end
    
    XBase = mean(XWin);
    XNorm = bsxfun(@times,XWin,1./XBase)-1;
    
    % POS
    Xs = squeeze(XNorm(:,2)-XNorm(:,3));
    Ys = squeeze(-2*XNorm(:,1)+XNorm(:,2)+XNorm(:,3));
    SWin = Xs + std(Xs)./std(Ys).*Ys;
    
    SWin = SWin - mean(SWin);
    %overlap, add
    if(i+WinL-1>length(S))%end of signal
        S(i:end)=S(i:end)+SWin(1:length(S(i:end)));
    else
        S(i:i+WinL-1)=S(i:i+WinL-1)+SWin;
    end
end

BVP_N=S-mean(S); %normalize
BVP=BVP_N;

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
    xlim([T(1) T(end)])%to account for differences caused by windows
end%endif plot

%% Remove Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
    rmpath([cd '\optional\']);
end

end%end function