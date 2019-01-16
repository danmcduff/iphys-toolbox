function [BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF)
% POS_WANG The Plane Orthogonal to Skin-Tone (POS) Method from: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491. DOI: 10.1109/TBME.2016.2609282
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
% Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

%% Parameters
SkinSegmentTF = false;

LPF = 0.7; %low cutoff frequency (Hz) - specified as 40 bpm (~0.667 Hz) in reference
HPF = 2.5; %high cutoff frequency (Hz) - specified as 240 bpm in reference

WinSec=1.6;%(based on refrence's 32 frame window with a 20 fps camera)

%% Add Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)
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

FramesToRead = floor(Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Read Video and Spatially Average:
T = zeros(FramesToRead,1);%initialize time vector
RGB = zeros(FramesToRead,3);%initialize color signal
FN = 0;
while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in
    %reference as a CSK detector from Henriques et al., 2012
    VidROI = VidFrame;
    
    if(SkinSegmentTF)%skin segmentation - originally specified in reference as an OC-SVM from Wang et al. 2015
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
    else
        RGB(FN,:) = sum(sum(VidROI,2)) ./ (size(VidROI,1)*size(VidROI,2));
    end
end

%% POS:
% Transform from: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017, May). Color-distortion filtering for remote photoplethysmography. In Automatic Face & Gesture Recognition (FG 2017), 2017 12th IEEE International Conference on (pp. 71-78). IEEE.
useFGTransform=false;
if useFGTransform
    RGBBase = mean(RGB);
    RGBNorm = bsxfun(@times,RGB,1./RGBBase)-1;
    FF = fft(RGBNorm);
    F = (0:size(RGBNorm,1)-1)*FS/size(RGBNorm,1);
    H = FF*[-1/sqrt(6);2/sqrt(6);-1/sqrt(6)];
    W = (H.*conj(H))./sum(FF.*conj(FF),2);
    FMask = (F >= LPF)&(F <= HPF);
    % FMask(length(FMask)/2+1:end)=FMask(length(FMask)/2:-1:1);
    FMask = FMask + fliplr(FMask);
    W=W.*FMask';%rectangular filter in frequency domain - not specified in original paper
    FF = FF.*repmat(W,[1,3]);
    RGBNorm=real(ifft(FF));
    RGBNorm = bsxfun(@times,RGBNorm+1,RGBBase);
    
    RGB=RGBNorm;
end

%lines and comments correspond to pseudo code algorithm on reference page 7       
N = size(RGB,1);%line 0 - RGB is of length N frames
H = zeros(1,N);%line 1 - initialize to zeros of length of video sequence
l = ceil(WinSec*FS);%line 1 - window length equivalent to reference: 32 samples of 20 fps camera (default 1.6s)
C = zeros(length(l),3);
for n = 1:N-1%line 2 - loop from first to last frame in video sequence
    %line 3 - spatial averaging was performed when video was read
    m = n - l + 1;%line 4 condition
    if(m > 0)%line 4
        Cn = ( RGB(m:n,:) ./ mean(RGB(m:n,:)) )';%line 5 - temporal normalization
        S = [0, 1, -1; -2, 1, 1] * Cn;%line 6 - projection
        h = S(1,:) + ((std(S(1,:)) / std(S(2,:))) * S(2,:));%line 7 - tuning
        H(m:n) = H(m:n) + (h - mean(h));%line 8 - overlap-adding
    end%line 9 - end if
end%line 10 - end for

BVP=H;
T=T(1:length(BVP));

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
        title('POS Method - ECG, PPG, iPPG Timeseries')
        
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
        title('POS Method - ECG, iPPG Timeseries')
        
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
