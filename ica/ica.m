function [W, Zhat] = ica(X,Nsources,Wprev)
%ICA  Perform independent component analysis
%  [W, ZHAT] = ICA(X) performs independent component analysis on data
%  observation matrix X. Matrix X is a transposed observation matrix, such that
%  each row of X represents an observed signal. This approach uses Cardoso's ICA
%  algorithm to estimate sources (ZHAT) and the de-mixing matrix W, an
%  approximation to A^{-1}, the original (unknown) mixing matrix. 

% Created by: G D Clifford 2004  gari AT mit DOT edu
% Last Modified: 5/7/06, documentation updated.

% Input argument checking
%------------------------
[a, b] = size(X);
if a > b
    fprintf('Warning - ICA works across the rows of the input data.\n');
    error('Please transpose input.');
end
%Nsources = a;

if Nsources > min([a b])
    Nsources = min([a b]);
    fprintf('Warning - number of soures cannot exceed number of observation channels \n')
    fprintf(' ... reducing to %i \n',Nsources)
end

%tic

if exist('Vprev','var'), [Winv, Zhat] = jade(X,Nsources,Wprev); else, [Winv, Zhat] = jade(X,Nsources); end
W = pinv(Winv);
%fprintf('algorithm timing ...  ')
%toc
