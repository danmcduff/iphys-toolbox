function [BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, FS, StartTime, Duration, ECGData, PPGData, PlotTF)
% POS_WANG The POS (Wang et al. 2016) Method.
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
SkinSegmentTF=false;

LPF = 0.70; %low cutoff frequency (Hz)
HPF = 4.0; %high cutoff frequency (Hz)
% NOTE DIFFERENT FROM AS IN THE PAPER FOR CONSISTENCY:
%LPF = 0.6667; %low cutoff frequency (Hz) - specified as 40 bpm.
%HPF = 4.0; %high cutoff frequency (Hz)

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
end

%% POS:
% Transform from: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017, May). Color-distortion filtering for remote photoplethysmography. In Automatic Face & Gesture Recognition (FG 2017), 2017 12th IEEE International Conference on (pp. 71-78). IEEE.
useFGTransform=0;
if useFGTransform
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
else
    RGBNorm = RGB;
end
        
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
