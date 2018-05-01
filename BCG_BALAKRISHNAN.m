function [BVP, PR, HR_ECG, PR_PPG, SNR] = BCG_BALAKRISHNAN(VideoFile, StartTime, Duration, BioSemiData, ECGMark, PPGMark, PlotTF)
% BCG_BALAKRISHNAN The BCG (Balakrishnan et al., 2012) Method Applied to the AFRL dataset.
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
bbox = [0,0,vidObj.width,vidObj.height];

minQ = 0.01;
points0 = detectMinEigenFeatures(rgb2gray(vidFrame),'ROI',bbox,'MinQuality',minQ);

% Remove eye region:
for i = size(points0.Location,1):-1:1
    if points0.Location(i,2)>= bbox(2)+round(boxH*0.2) && points0.Location(i,2)<= bbox(2)+round(boxH*0.55)
        points0(i) = [];
    end
end

if length(points0) < 5
    HR = 0; HR0 = 0; mySNR = 0;
    return
end

%% Intitalize:
d = [];
d = [d points0.Location(:,2)];
tracker = vision.PointTracker;
initialize(tracker,points0.Location,vidFrame);

%% Initialize Time Vector:
t = zeros(framesToRead,1);

%% Read Video and Spatially Average:
fn=0;
while hasFrame(vidObj) && (vidObj.CurrentTime <= StartTime+Duration)
    fn = fn+1;
    t(fn) = vidObj.CurrentTime;
    vidFrame = readFrame(vidObj);
    [points, validity] = step(tracker,vidFrame);
    d = [d points(:,2)];
end

%% Butterworth 3rd order filter:
lpf = 0.70; %low cutoff frequency (Hz)
hpf = 2.50; %high cutoff frequency (Hz)
[b,a] = butter(3,[lpf/nyquist hpf/nyquist]);
d_filt = filtfilt(b,a,double(d'));

d_filt2 = bsxfun(@minus, d_filt, d_filt(1,:));
dL2 = sqrt(sum(d_filt2.^2,2));
dMask = dL2<=quantile(dL2,.75);
d_filt3=d_filt2(dMask,:);
[coeff, ~, latent] = pca(d_filt3);
score = bsxfun(@minus, d_filt, mean(d_filt))/coeff';

[pxx2,f2] = plomb(score(:,1:5),t);

fmask2 = (f2 >= 0.7)&(f2 <= 2.5);
frange2 = f2(fmask2);
R = zeros(1,5);
for i = 1:5
    [maxP,idxP] = max(pxx2(:,i));
    R(i) = (maxP+pxx2(idxP*2,i))/sum(pxx2(:,i));
end
[~,numPC] = max(R);
PR = frange2(argmax(pxx2(fmask2,numPC),1))*60;

BVP_f= score(:,numPC);

BVP_n=BVP_f-mean(BVP_f); 
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

