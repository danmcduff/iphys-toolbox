function [BVP, PR, HR_ECG, PR_PPG, SNR] = POS_WANG(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% POS_WANG The POS (Wang et al. 2016) Method Applied to the AFRL dataset.
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
%       HR_ECG                  = Gold startdard Heart Rate measured from the ECG timeseries for the window.
%       PR_PPG                  = Pulse Rate measured from the PPG timeseries for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, April 2018

%% Plot Control
if(PlotTF)
    PlotPRPSD = true;
    PlotSNR = true;
else
    PlotPRPSD = false;
    PlotSNR = false;
end

%% Load Video:
vidObj = VideoReader(VideoFile);
vidObj.CurrentTime = StartTime;

fs = 120; %true frame rate
framesToRead=ceil(Duration*vidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Initialize Face Detector (UNCOMMENT IF USING A FACE DETECTOR):
% faceDetector = vision.CascadeObjectDetector();
% d_raw=[];
% bbox = step(faceDetector, vidFrame);
% [~,Idx]=max(bbox(:,3));
% bbox=bbox(Idx,:);

%% Initialize Time Vector:
t = zeros(framesToRead,1);

%% Read Video and Spatially Average:
colors=zeros(framesToRead,3);
fn=0;
while hasFrame(vidObj) && (vidObj.CurrentTime <= StartTime+Duration)
    fn = fn+1;
    t(fn) = vidObj.CurrentTime;
    vidFrame = readFrame(vidObj);
    vidROI = vidFrame;
    
    vidROI = vidFrame;
    YCBCR = rgb2ycbcr(vidROI);
    Yth = YCBCR(:,:,1)>80;
    CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
    CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
    ROISkin = vidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
    d_raw = [d_raw sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2)]; 
end
colors = sum(sum(vidROI));

%% POS:
d = colors;
N = size(d,1);
winL = 32/30*120; % /20*120
Nwin = N-winL+1;
S = zeros(Nwin,1);
for i = 1:Nwin
    Xwin = d(i:i+winL-1,:);
    Xbase = mean(Xwin);
    Xnorm = bsxfun(@times,Xwin,1./Xbase)-1;
    
    % CHROM
    Xs = squeeze(3*Xnorm(:,1)-2*Xnorm(:,2));
    Ys = squeeze(1.5*Xnorm(:,1)+Xnorm(:,2)-1.5*Xnorm(:,3));
    Swin = Xs - std(Xs)./std(Ys).*Ys;

    Swin = Swin - mean(Swin);
    S(i)=Swin(round(winL/2));
end

BVP_i = S; 
BVP_f = BVP_i;

BVP_n=BVP_f-mean(BVP_f); %normalize
BVP=BVP_n;

% Estimate Pulse Rate from periodogram
PR = prpsd(BVP,fs,lpf*60,hpf*60,PlotPRPSD);

%% Groundtruth HR:
ECG=load(ECGMark,'KeepBeatData');
%KeepBeatData Format (ECG)- Beat Index, Inter-Beat Interval (ms), Zero-Referenced Sample Index of R-wave Maximum, Zero-Referenced Time of R-Wave (s), Amplitude of R-Wave in Filtered Signal, Zero-Referenced Sample Index of R-wave Maximum, Zero-Referenced Time of S-Wave (s), Amplitude of S-Wave in Filtered Signal, User Added Data (0-automatic, 1-user added/ corrected)
ECGMask = (ECG.KeepBeatData(:,4)>=StartTime)&(ECG.KeepBeatData(:,4)<=StartTime+Duration);
HR_ECG = 1/mean(diff(ECG.KeepBeatData(ECGMask,4)))*60;

PPG=load(PPGMark,'KeepBeatData');
%KeepBeatData Format (PPG)- Pulse Index, Inter-Pulse Interval (ms), Zero-Referenced Sample Index of Systolic Onset, Zero-Referenced Time of Systolic Onset (s), Amplitude of Pulse/ Systolic Onset in Filtered Signal, Empty, Empty, Empty, User Added Data (0-automatic, 1-user added/ corrected)
PPGMask = (PPG.KeepBeatData(:,4)>=StartTime)&(PPG.KeepBeatData(:,4)<=StartTime+Duration);
PR_PPG = 1/mean(diff(PPG.KeepBeatData(PPGMask,4)))*60;

%% SNR
SNR = bvpsnr(BVP,fs,HR_ECG,PlotSNR);

%% Optionally Plot Timeseries
if(PlotTF)
    %{
    % Plot iPPG BVP
    figure, plot(t,BVP),xlabel('Time (s)'),ylabel('Amplitude (a.u.)')
    %}
    
    %Plot ECG, PPG, iPPG timeseries
    load(BioSemiData,'FiltNECG', 'FiltNPPG', 'ResampleTimeSeries')%load filtered and resampled ECG, PPG timeseries data
    figure
    subplot(3,1,1)
    plot(ResampleTimeSeries,FiltNECG)
    hold on
    plot(ECG.KeepBeatData(:,4),ECG.KeepBeatData(:,5),'*')
    xlim([StartTime StartTime+Duration])
    ylabel('ECG (a.u.)')
    title('ECG, PPG, iPPG Timeseries')
    
    subplot(3,1,2)
    plot(ResampleTimeSeries,FiltNPPG)
    hold on
    plot(PPG.KeepBeatData(:,4),-1*PPG.KeepBeatData(:,5),'*')
    xlim([StartTime StartTime+Duration])
    ylabel('PPG (a.u.)')
    
    subplot(3,1,3)
    plot(t,BVP)
    hold on
    xlim([StartTime StartTime+Duration])
    ylabel('iPPG (a.u.)')
    
    xlabel('Time (s)')
end
end

