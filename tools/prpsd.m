function [PR] = prpsd(BVP, fs, llPR, ulPR, PlotTF)
%PRPSD Estimates a pulse rate from a BVP signal
%   Inputs:
%       BVP                     = A BVP timeseries.
%       fs                      = The sample rate of the BVP time series (Hz/fps).
%       llPR                    = The lower limit for pulse rate (bpm).
%       ulPR                    = The upper limit for pulse rate (bpm).
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       PR                      = The estimated PR in BPM.
% Daniel McDuff, Ethan Blackford, Justin Estepp, April 2018
%%
nyquist = fs/2;
fResBPM = 0.5; %resolution (bpm) of bins in power spectrum used to determine PR and SNR

N = (60*2*nyquist)/fResBPM;

%% Construct Periodogram
[pxx,f] = periodogram(BVP,hamming(length(BVP)),N,fs);
fMask = (f >= (llPR/60))&(f <= (ulPR/60));

%% Calculate predicted HR:
fRange = f(fMask);
pRange = pxx(fMask);
maxInd = argmax(pxx(fMask),1);
PR_f = fRange(maxInd);
PR = PR_f*60;

%% Optionally Plot the PSD and peak frequency
if(PlotTF)
    %{
    % Plot PSD and peak frequency
    figure
    plot(f,pxx)
    hold on
    plot(PR_f,pRange(maxInd),'*r')
    text(PR_f,pRange(maxInd),['   ' num2str(PR_f,'%3.2f') ' Hz; ' num2str(PR,'%4.1f'),' bpm'])
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)')
    xlim([0 4])
    title('Power Spectrum and Peak Frequency')
    %}
    
    % Plot PSD (in dB) and peak frequency
    figure
    plot(f,pow2db(pxx))
    hold on
    plot(PR_f,pow2db(pRange(maxInd)),'*r')
    text(PR_f,pow2db(pRange(maxInd)),['   ' num2str(PR_f,'%3.2f') ' Hz; ' num2str(PR,'%4.1f'),' bpm'])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    xlim([0 4])
    title('Power Spectrum and Peak Frequency')
end
end

