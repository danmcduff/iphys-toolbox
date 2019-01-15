function [W, Zhat] = ica(X,Nsources,Wprev)
%ICA  Perform independent component analysis
%  [W, ZHAT] = ICA(X) performs independent component analysis on data
%  observation matrix X of physiological observations.
%
%   Inputs:
%       X               = Observation matrix: rows should be observations, cols samples.
%       Nsources        = Number of source signals (optional).
%       Wprev           = Timepoint at which to start process (default = 0 seconds).
%
%   Outputs:
%       W               = Demixing matrix.
%       Zhat            = Source signal matrix.
%
%  This approach uses Cardoso's ICA
%  algorithm to estimate sources (ZHAT) and the de-mixing matrix W, an
%  approximation to A^{-1}, the original (unknown) mixing matrix. 
%
% Daniel McDuff, Ethan Blackford, Justin Estepp, December 2018
% Based on an implementation by: G D Clifford (2004) gari AT mit DOT edu
% Licensed under the RAIL AI License.

[nRows, nCols] = size(X);
if nRows > nCols
    fprintf('Warning - The number of rows is cannot be greater than the number of columns.\n');
    error('Please transpose input.');
end

if Nsources > min([nRows nCols])
    Nsources = min([nRows nCols]);
    fprintf('Warning - The number of soures cannot exceed number of observation channels. \n')
    fprintf('The number of sources will be reduced to the number of observation channels (%i) \n',Nsources)
end

if exist('Vprev','var')
    [Winv, Zhat] = jade(X,Nsources,Wprev); 
else
    [Winv, Zhat] = jade(X,Nsources); 
end
W = pinv(Winv);
