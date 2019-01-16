function [PR] = prpsd(BVP, FS, LL_PR, UL_PR, PlotTF)
%PRPSD Estimates a pulse rate from a BVP signal
%   Inputs:
%       BVP                     = A BVP timeseries.
%       FS                      = The sample rate of the BVP time series (Hz/fps).
%       LL_PR                    = The lower limit for pulse rate (bpm).
%       UL_PR                    = The upper limit for pulse rate (bpm).
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       PR                      = The estimated PR in BPM.
%
% Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

%%
Nyquist = FS/2;
FResBPM = 0.5; %resolution (bpm) of bins in power spectrum used to determine PR and SNR

N = (60*2*Nyquist)/FResBPM;

%% Construct Periodogram
[Pxx,F] = periodogram(BVP,hamming(length(BVP)),N,FS);
FMask = (F >= (LL_PR/60))&(F <= (UL_PR/60));

%% Calculate predicted HR:
FRange = F(FMask);
PRange = Pxx(FMask);
MaxInd = argmax(Pxx(FMask),1);
PR_F = FRange(MaxInd);
PR = PR_F*60;

%% Optionally Plot the PSD and peak frequency
if(PlotTF)
    %{
    % Plot PSD and peak frequency
    figure
    plot(F,Pxx)
    hold on
    plot(PR_F,pRange(MaxInd),'*r')
    text(PR_F,pRange(MaxInd),['   ' num2str(PR_F,'%3.2f') ' Hz; ' num2str(PR,'%4.1f'),' bpm'])
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)')
    xlim([0 4.5])
    title('Power Spectrum and Peak Frequency')
    %}
    
    % Plot PSD (in dB) and peak frequency
    figure
    plot(F,pow2db(Pxx))
    hold on
    plot(PR_F,pow2db(PRange(MaxInd)),'*r')
    text(PR_F,pow2db(PRange(MaxInd)),['   ' num2str(PR_F,'%3.2f') ' Hz; ' num2str(PR,'%4.1f'),' bpm'])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    xlim([0 4.5])
    title('Power Spectrum and Peak Frequency')
end
end

