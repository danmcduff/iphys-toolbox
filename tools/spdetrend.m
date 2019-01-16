function [y] = spdetrend(x,lambda)
%SPDETREND Detrends a signal with the smoothness priors approach implemented 
% based on M. P. Tarvainen, TBME, 2002
%
%   Inputs:
%       x                       = Signal to detrend.
%       lambda                  = lambda value to use for detrending.
%
%   Outputs:
%       y                       = The detrended signal.
%
% Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

T=length(x);
I=speye(T);
D2=spdiags(ones(T-2,1)*[1 -2 1],[0:2],T-2,T);
sr=double(x);
y=(I-inv(I+lambda^2*(D2'*D2)))*sr;
end

