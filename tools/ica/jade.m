function [A,S]=jade(X,m,Wprev)
% ICA implementation using JADE (used by the ICA function).
%
%   Inputs:
%       X               = Each column of X is a sample from the n sensors.
%       m               = Number of source signals (optional).
%       Wprev           = Timepoint at which to start process (default = 0 seconds).
%
%   Outputs:
%       A               = An n x m estimate of an unknown matrix.
%       S               = An m x T naive estimate of the source signals.
%
% Implementation by: JF Cardoso. (Version 1.5)
% 
% Please Cite: 
% Reference:
%  @article{CS_iee_94,
%   author = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
%   journal = "IEE Proceedings-F",
%   title = "Blind beamforming for non {G}aussian signals",
%   number = "6",
%   volume = "140",
%   month = dec,
%   pages = {362-370},
%   year = "1993"}

% Determine the no of rows and colums of the observation matrix.
[n,T]	= size(X);
if nargin==1, m=n ; end;

nem     = m;                    % The nuumber of eigen-matrices to be diagonalized.
seuil	= 1/sqrt(T)/100;    % A statistical threshold for stopping the joint diag.

%% Whiten the Matrix:
% Assuming white noise:
if m<n, 
 	[U,D]       = eig((X*X')/T); 
	[pu,k]   = sort(diag(D));
 	ibl         = sqrt(pu(n-m+1:n)-mean(pu(1:n-m)));
 	bl          = ones(m,1) ./ ibl ;
 	W           = diag(bl)*U(1:n,k(n-m+1:n))';
 	IW          = U(1:n,k(n-m+1:n))*diag(ibl);
% Assuming no noise:
else 
 	IW          = sqrtm((X*X')/T);
 	W           = inv(IW);
end;
Y	= W*X;

%% Cumulant Estimation:
R       = (Y*Y' )/T ;
C       = (Y*Y.')/T ;
Q       = zeros(m*m*m*m,1) ;
index	= 1;

for lx = 1:m ; Yl 	= Y(lx,:);
    for kx = 1:m ; Ykl 	= Yl.*conj(Y(kx,:));
        for jx = 1:m ; Yjkl	= Ykl.*conj(Y(jx,:));
            for ix = 1:m ; 
                Q(index) = ...
                (Yjkl * Y(ix,:).')/T -  R(ix,jx)*R(lx,kx) -  R(ix,kx)*R(lx,jx) -  C(ix,lx)*conj(C(jx,kx))  ;
                index	= index + 1;
            end
        end
    end
end

%% Compute and Reshape the significant Eigen Matrices:
[U,D]	= eig(reshape(Q,m*m,m*m)); 
[la,K]	= sort(abs(diag(D)));

M	= zeros(m,nem*m);	
Z	= zeros(m)	;
h	= m*m;
for u=1:m:nem*m, 
	Z(:)            = U(:,K(h));
	M(:,u:u+m-1)	= la(h)*Z;
	h               = h-1; 
end;

%% Approximate the Diagonalization of the Eigen Matrices:
B       = [ 1 0 0 ; 0 1 1 ; 0 -i i ] ;
Bt      = B' ;

encore	= 1;
if exist('Wprev','var'), V=inv(Wprev); else, V	= eye(m); end

%% Main Loop:
while encore, encore=0;
    for p = 1:m-1,
        for q = p+1:m,
            
            Ip = p:m:nem*m ;
            Iq = q:m:nem*m ;
            
            % Computing the Givens Angles:
            g	= [ M(p,Ip)-M(q,Iq)  ; M(p,Iq) ; M(q,Ip) ] ;
            [vcp,D] = eig(real(B*(g*g')*Bt));
            [la, K]	= sort(diag(D));
            angles	= vcp(:,K(3));
            if angles(1)<0 , angles= -angles ; end ;
            c	= sqrt(0.5+angles(1)/2);
            s	= 0.5*(angles(2)-j*angles(3))/c;
            
            % Update the matrices M and V by a Givens Rotation:
            if abs(s)>seuil,
                encore          = 1 ;
                pair            = [p;q] ;
                G               = [ c -conj(s) ; s c ] ;
                V(:,pair)       = V(:,pair)*G ;
                M(pair,:)       = G' * M(pair,:) ;
                M(:,[Ip Iq]) 	= [ c*M(:,Ip)+s*M(:,Iq) -conj(s)*M(:,Ip)+c*M(:,Iq) ] ;
            end
        end
    end
end

%% Estimation of the mixing matrix and signal separation
A	= IW*V;
S	= V'*Y ;

return ;

