function [SNR] = bvpsnr(BVP, FS, HR, PlotTF)
%BVPSNR Estimates the signal-to-noise ratio of the blood volume pulse
% signal. Adapted from the method by G. de Haan, TBME, 2013.
% SNR calculated as the ratio (in dB) of power contained within +/- 0.1 Hz
% of the reference heart rate frequency and +/- 0.2 of its first
% harmonic and sum of all other power between 0.5 and 4 Hz.
% Adapted from the method by G. de Haan, TBME, 2013
%
%   Inputs:
%       BVP                     = A BVP timeseries.
%       FS                      = The sample rate of the BVP time series (Hz/fps).
%       HR                      = The reference heart rate (Hz/fps).
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio.
% Daniel McDuff, Ethan Blackford, Justin Estepp, June 2018

HR_F=HR/60;

NyquistF = FS/2;
FResBPM = 0.5; %resolution (bpm) of bins in power spectrum used to determine PR and SNR
N = (60*2*NyquistF)/FResBPM; %number of bins in power spectrum

%% Construct Periodogram
[Pxx,F] = periodogram(BVP,hamming(length(BVP)),N,FS);
GTMask1 = (F >= HR_F-0.1)&(F <= HR_F+0.1);
GTMask2 = (F >= HR_F*2-0.2)&(F <= HR_F*2+0.2);
SPower = sum(Pxx(GTMask1|GTMask2));
FMask2 = (F >= 0.5)&(F <= 4);
AllPower = sum(Pxx(FMask2));
SNR = pow2db(SPower/(AllPower-SPower));

%% Optionally plot the power spectrum and regions used to calculate the SNR
if(PlotTF)
    % Plot spower spectrum and SNR regions
    figure
    plot(F,pow2db(Pxx))
    title('Power Spectrum and SNR Regions')
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    xlim([0.5 4])
    ylimreg=ylim;
    hold on
    
    % HR peak
    %line([HR_F,HR_F],ylim,'Color','magenta','LineStyle','-.')
    
    % HR peak region
    line([HR_F-0.1,HR_F-0.1],ylim,'Color','red','LineStyle','--')
    line([HR_F+0.1,HR_F+0.1],ylim,'Color','red','LineStyle','--')
    
    % First harmonic
    line([HR_F*2-0.2,HR_F*2-0.2],ylim,'Color','red','LineStyle','--')
    line([HR_F*2+0.2,HR_F*2+0.2],ylim,'Color','red','LineStyle','--')
    
    % Overall power region
    line([0.5,0.5],ylim,'Color','black','LineStyle','-')
    line([4,4],ylim,'Color','black','LineStyle','-')
    
    % Adjust axes
    xlim([0 4.5])
    ylim(ylimreg);
end

end

