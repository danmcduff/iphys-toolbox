function [SNR] = bvpsnr(BVP, fs, HR, PlotTF)
%BVPSNR Estimates the signal-to-noise ratio of the blood volume pulse
% signal. Adapted from the method by G. de Haan, TBME, 2013.
% SNR calculated as the ratio (in dB) of power contained within +/- 0.1 Hz
% of the reference heart rate frequency and +/- 0.2 of its first
% harmonic and sum of all other power between 0.5 and 4 Hz.
% Adapted from the method by G. de Haan, TBME, 2013
%
%   Inputs:
%       BVP                     = A BVP timeseries.
%       fs                      = The sample rate of the BVP time series (Hz/fps).
%       HR                      = The reference heart rate (Hz/fps).
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio.
% Daniel McDuff, Ethan Blackford, Justin Estepp, April 2018

HRf=HR/60;

nyquist = fs/2;
fResBPM = 0.5; %resolution (bpm) of bins in power spectrum used to determine PR and SNR
N = (60*2*nyquist)/fResBPM; %number of bins in power spectrum

%% Construct Periodogram
[pxx,f] = periodogram(BVP,hamming(length(BVP)),N,fs);
gtMask1 = (f >= HRf-0.1)&(f <= HRf+0.1);
gtMask2 = (f >= HRf*2-0.2)&(f <= HRf*2+0.2);
sPower = sum(pxx(gtMask1|gtMask2));
fMask2 = (f >= 0.5)&(f <= 4);
allPower = sum(pxx(fMask2));
SNR = pow2db(sPower/(allPower-sPower));

%% Optionally plot the power spectrum and regions used to calculate the SNR
if(PlotTF)
    % Plot spower spectrum and SNR regions
    figure
    plot(f,pow2db(pxx))
    title('Power Spectrum and SNR Regions')
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    xlim([0.5 4])
    ylimreg=ylim;
    hold on
    
    % HR peak
    %line([HRf,HRf],ylim,'Color','magenta','LineStyle','-.')
    
    % HR peak region
    line([HRf-0.1,HRf-0.1],ylim,'Color','red','LineStyle','--')
    line([HRf+0.1,HRf+0.1],ylim,'Color','red','LineStyle','--')
    
    % First harmonic
    line([HRf*2-0.2,HRf*2-0.2],ylim,'Color','red','LineStyle','--')
    line([HRf*2+0.2,HRf*2+0.2],ylim,'Color','red','LineStyle','--')
    
    % Overall power region
    line([0.5,0.5],ylim,'Color','black','LineStyle','-')
    line([4,4],ylim,'Color','black','LineStyle','-')
    
    % Adjust axes
    xlim([0 4.5])
    ylim(ylimreg);
end

end

