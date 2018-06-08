function [A,S]=jade(X,m,Vprev)
% Source separation of complex signals with JADE.
% Jade performs `Source Separation' in the following sense:
%   X is an n x T data matrix assumed modelled as X = A S + N where
% 
% o A is an unknown n x m matrix with full rank.
% o S is a m x T data matrix (source signals) with the properties
%    	a) for each t, the components of S(:,t) are statistically
%    	   independent
% 	b) for each p, the S(p,:) is the realization of a zero-mean
% 	   `source signal'.
% 	c) At most one of these processes has a vanishing 4th-order
% 	   cumulant.
% o  N is a n x T matrix. It is a realization of a spatially white
%    Gaussian noise, i.e. Cov(X) = sigma*eye(n) with unknown variance
%    sigma.  This is probably better than no modeling at all...
%
% Jade performs source separation via a 
% Joint Approximate Diagonalization of Eigen-matrices.  
%
% THIS VERSION ASSUMES ZERO-MEAN SIGNALS
%
% Input :
%   * X: Each column of X is a sample from the n sensors
%   * m: m is an optional argument for the number of sources.
%     If ommited, JADE assumes as many sources as sensors.
%
% Output :
%    * A is an n x m estimate of the mixing matrix
%    * S is an m x T naive (ie pinv(A)*X)  estimate of the source signals
%
%
% Version 1.5.  Copyright: JF Cardoso.  
%
% See notes, references and revision history at the bottom of this file



[n,T]	= size(X);

%%  source detection not implemented yet !
if nargin==1, m=n ; end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A few parameters that could be adjusted
nem	= m;		% number of eigen-matrices to be diagonalized
seuil	= 1/sqrt(T)/100;% a statistical threshold for stopping joint diag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% whitening
%
if m<n, %assumes white noise
 	[U,D] 	= eig((X*X')/T); 
	[puiss,k]=sort(diag(D));
 	ibl 	= sqrt(puiss(n-m+1:n)-mean(puiss(1:n-m)));
 	bl 	= ones(m,1) ./ ibl ;
 	W	= diag(bl)*U(1:n,k(n-m+1:n))';
 	IW 	= U(1:n,k(n-m+1:n))*diag(ibl);
else    %assumes no noise
 	IW 	= sqrtm((X*X')/T);
 	W	= inv(IW);
end;
Y	= W*X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cumulant estimation


R	= (Y*Y' )/T ;
C	= (Y*Y.')/T ;

Yl	= zeros(1,T);
Ykl	= zeros(1,T);
Yjkl	= zeros(1,T);

Q	= zeros(m*m*m*m,1) ;
index	= 1;

for lx = 1:m ; Yl 	= Y(lx,:);
for kx = 1:m ; Ykl 	= Yl.*conj(Y(kx,:));
for jx = 1:m ; Yjkl	= Ykl.*conj(Y(jx,:));
for ix = 1:m ; 
	Q(index) = ...
	(Yjkl * Y(ix,:).')/T -  R(ix,jx)*R(lx,kx) -  R(ix,kx)*R(lx,jx) -  C(ix,lx)*conj(C(jx,kx))  ;
	index	= index + 1 ;
end ;
end ;
end ;
end

%% If you prefer to use more memory and less CPU, you may prefer this
%% code (due to J. Galy of ENSICA) for the estimation the cumulants
%ones_m = ones(m,1) ; 
%T1 	= kron(ones_m,Y); 
%T2 	= kron(Y,ones_m);  
%TT 	= (T1.* conj(T2)) ;
%TS 	= (T1 * T2.')/T ;
%R 	= (Y*Y')/T  ;
%Q	= (TT*TT')/T - kron(R,ones(m)).*kron(ones(m),conj(R)) - R(:)*R(:)' - TS.*TS' ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%computation and reshaping of the significant eigen matrices

[U,D]	= eig(reshape(Q,m*m,m*m)); 
[la,K]	= sort(abs(diag(D)));

%% reshaping the most (there are `nem' of them) significant eigenmatrice
M	= zeros(m,nem*m);	% array to hold the significant eigen-matrices
Z	= zeros(m)	; % buffer
h	= m*m;
for u=1:m:nem*m, 
	Z(:) 		= U(:,K(h));
	M(:,u:u+m-1)	= la(h)*Z;
	h		= h-1; 
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% joint approximate diagonalization of the eigen-matrices


%% Better declare the variables used in the loop :
B 	= [ 1 0 0 ; 0 1 1 ; 0 -i i ] ;
Bt	= B' ;
Ip	= zeros(1,nem) ;
Iq	= zeros(1,nem) ;
g	= zeros(3,nem) ;
G	= zeros(2,2) ;
vcp	= zeros(3,3);
D	= zeros(3,3);
la	= zeros(3,1);
K	= zeros(3,3);
angles	= zeros(3,1);
pair	= zeros(1,2);
c	= 0 ;
s	= 0 ;


%init;
encore	= 1;
if exist('Vprev','var'), V=inv(Vprev); else, V	= eye(m); end

% Main loop
while encore, encore=0;
 for p=1:m-1,
  for q=p+1:m,

 	Ip = p:m:nem*m ;
	Iq = q:m:nem*m ;

	% Computing the Givens angles
 	g	= [ M(p,Ip)-M(q,Iq)  ; M(p,Iq) ; M(q,Ip) ] ; 
 	[vcp,D] = eig(real(B*(g*g')*Bt));
	[la, K]	= sort(diag(D));
 	angles	= vcp(:,K(3));
	if angles(1)<0 , angles= -angles ; end ;
 	c	= sqrt(0.5+angles(1)/2);
 	s	= 0.5*(angles(2)-j*angles(3))/c; 

 	if abs(s)>seuil, %%% updates matrices M and V by a Givens rotation
	 	encore 		= 1 ;
		pair 		= [p;q] ;
 		G 		= [ c -conj(s) ; s c ] ;
		V(:,pair) 	= V(:,pair)*G ;
	 	M(pair,:)	= G' * M(pair,:) ;
		M(:,[Ip Iq]) 	= [ c*M(:,Ip)+s*M(:,Iq) -conj(s)*M(:,Ip)+c*M(:,Iq) ] ;
 	end%% if
  end%% q loop
 end%% p loop
end%% while

%%%estimation of the mixing matrix and signal separation
A	= IW*V;
S	= V'*Y ;

return ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Note 1: This version does *not* assume circularly distributed
% signals as 1.1 did.  The difference only entails more computations
% in estimating the cumulants
%
%
% Note 2: This code tries to minimize the work load by jointly
% diagonalizing only the m most significant eigenmatrices of the
% cumulant tensor.  When the model holds, this avoids the
% diagonalization of m^2 matrices.  However, when the model does not
% hold, there is in general more than m significant eigen-matrices.
% In this case, this code still `works' but is no longer equivalent to
% the minimization of a well defined contrast function: this would
% require the diagonalization of *all* the eigen-matrices.  We note
% (see the companion paper) that diagonalizing **all** the
% eigen-matrices is strictly equivalent to diagonalize all the
% `parallel cumulants slices'.  In other words, when the model does
% not hold, it could be a good idea to diagonalize all the parallel
% cumulant slices.  The joint diagonalization will require about m
% times more operations, but on the other hand, computation of the
% eigen-matrices is avoided.  Such an approach makes sense when
% dealing with a relatively small number of sources (say smaller than
% 10).
%
%
% Revision history
%-----------------
%
% Version 1.5 (Nov. 2, 97) : 
% o Added the option kindly provided by Jerome Galy
%   (galy@dirac.ensica.fr) to compute the sample cumulant tensor.
%   This option uses more memory but is faster (a similar piece of
%   code was also passed to me by Sandip Bose).
% o Suppressed the useles variable `oui'.
% o Changed (angles=sign(angles(1))*angles) to (if angles(1)<0 ,
%   angles= -angles ; end ;) as suggested by Iain Collings
%   <i.collings@ee.mu.OZ.AU>.  This is safer (with probability 0 in
%   the case of sample statistics)
% o Cosmetic rewriting of the doc.  Fixed some typos and added new
%   ones.
%
% Version 1.4 (Oct. 9, 97) : Changed the code for estimating
% cumulants. The new version loops thru the sensor indices rather than
% looping thru the time index.  This is much faster for large sample
% sizes.  Also did some clean up.  One can now change the number of
% eigen-matrices to be jointly diagonalized by just changing the
% variable `nem'.  It is still hard coded below to be equal to the
% number of sources.  This is more economical and OK when the model
% holds but not appropriate when the model does not hold (in which
% case, the algorithm is no longer asymptotically equivalent to
% minimizing a contrast function, unless nem is the square of the
% number of sources.)
%
% Version 1.3 (Oct. 6, 97) : Added various Matalb tricks to speed up
% things a bit.  This is not very rewarding though, because the main
% computational burden still is in estimating the 4th-order moments.
% 
% Version 1.2 (Mar., Apr., Sept. 97) : Corrected some mistakes **in
% the comments !!**, Added note 2 `When the model does not hold' and
% the EUSIPCO reference.
%
% Version 1.1 (Feb. 94): Creation
%
%-------------------------------------------------------------------
%
% Contact JF Cardoso for any comment bug report,and UPDATED VERSIONS.
% email : cardoso@sig.enst.fr 
% or check the WEB page http://sig.enst.fr/~cardoso/stuff.html 
%
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
%
%
%  Some analytical insights into the asymptotic performance of JADE are in
% @inproceedings{PerfEusipco94,
%  HTML 	= "ftp://sig.enst.fr/pub/jfc/Papers/eusipco94_perf.ps.gz",
%  author       = "Jean-Fran\c{c}ois Cardoso",
%  address      = {Edinburgh},
%  booktitle    = "{Proc. EUSIPCO}",
%  month 	= sep,
%  pages 	= "776--779",
%  title 	= "On the performance of orthogonal source separation algorithms",
%  year 	= 1994}
%_________________________________________________________________________
% jade.m ends here
