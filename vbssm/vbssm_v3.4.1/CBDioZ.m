% Calculates the mean and variance of each element in CB+D.
% Also returns the significances of each entry.
%
%  [CBDmean,CBDvar,CBDZ,CBDZsig] = CBDioZ(net);
%
% MJB GNCU 02/07/03
% updated 19 Oct 05

function [CBDmean,CBDvar,CBDZ,CBDZsig] = CBDioZ(net,dispopt,sds);

if nargin<2, dispopt=0; end
if nargin<3, sds=3; end

% Preliminaries
[p k] = size(net.exp.C);
[pinp] = size(net.exp.D,2);
GC = net.param.MD;

% More preliminaries
Erhoinv = net.param.qrhob/(net.param.qrhoa-1);
SigC = net.param.SigC;
  SigD = net.param.SigD;
  SigCh = SigC+SigC*GC*SigD*GC'*SigC;
  C = net.exp.C;
  D = net.exp.D;
SigB = net.param.SigB;
  B = net.exp.B;
SC = net.param.SC; 

% Even more preliminaries
tmp = (C*B).^2;
tmp1 = SigC*SC;
tmp2 = SigC*GC*D';
tmp3 = SigC*GC*SigD;

% Now labouriously calculate variance on each element under the variational
% approximation.  See me to 

for ii = 1:p;
  if dispopt, fprintf('.'); end
  for jj = 1:pinp;
    
    var1tmp = zeros(k,k);
    for kk = 1:k;
      for ll = 1:k;
	var1tmp(kk,ll) = ( Erhoinv(ii)*SigCh(kk,ll)+C(ii,kk)*C(ii,ll) ...
				 )*( SigB(jj,jj)*(kk==ll) + B(kk,jj)*B(ll,jj) );
      end
    end
    var1(ii,jj) = sum(sum(var1tmp,1),2) - tmp(ii,jj);

    for kk = 1:k;
      zeta(ii,kk,jj) = tmp1(kk,ii)*D(ii,jj)-tmp2(kk,ii)*D(ii,jj)-tmp3(kk,jj)* ...
	  Erhoinv(ii);
      
      var3tmp(ii,kk,jj) = 2*( zeta(ii,kk,jj)-D(ii,jj)*C(ii,kk) ) * B(kk,jj);
    end % kk
  
  end % jj
end % ii

%var1
%done

%var2
var2 = zeros(p,pinp);
for ii = 1:p;
  for jj = 1:pinp;
    var2(ii,jj) = Erhoinv(ii)*SigD(jj,jj);
  end
end

%var3
var3 = squeeze( sum( var3tmp ,2) ); % i,j now

CBDvar = var1+var2+var3;  
CBDmean = C*B+D; % note these are means, and C and B are independent in the
                % posterior under the variational approximation.

CBDZ = CBDmean./sqrt(CBDvar);

if dispopt, fprintf('\n'); end

CBDZsig = CBDZ.*(abs(CBDZ)>sds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END





return

% BELOW is a previous attempt which gives the same result but at a slower
% pace!

% Calculates the variance of each element in CB+D.

% preliminaries
[p k] = size(net.exp.C);
[pinp] = size(net.exp.D,2);
GC = net.param.MD;

Erhoinv = net.param.qrhob/(net.param.qrhoa-1);
SigC = net.param.SigC;
  SigD = net.param.SigD;
  SigCh = SigC+SigC*GC*SigD*GC'*SigC;
  C = net.exp.C;
  D = net.exp.D;
SigB = net.param.SigB;
  B = net.exp.B;
SC = net.param.SC; 

tmp1 = SigC*SC;
tmp2 = SigC*GC*D';
tmp3 = SigC*GC*SigD;

var1tmp = zeros(p,pinp,k,k);
zeta = zeros(p,k,pinp);
var3tmp = zeros(p,k,pinp);
for ii = 1:p;
  if dispopt, fprintf('.'); end
  for jj = 1:pinp;
    for kk = 1:k;

      for ll = 1:k;
	var1tmp(ii,jj,kk,ll) = Erhoinv(ii)*SigCh(ll,kk)*SigB(jj,jj)*(ll==kk) ...
	    + Erhoinv(ii)*SigCh(ll,kk)*B(ll,jj)*B(kk,jj) ...
	    + C(ii,ll)*C(ii,kk)*SigB(jj,jj)*(ll==kk);
      end	

      zeta(ii,kk,jj) = tmp1(kk,ii)*D(ii,jj)-tmp2(kk,ii)*D(ii,jj)-tmp3(kk,jj)* ...
	  Erhoinv(ii);
      
      var3tmp(ii,kk,jj) = 2*( zeta(ii,kk,jj)-D(ii,jj)*C(ii,kk) ) * B(kk, jj);
    end % kk
  end % jj
end % ii

%var1
var1 = sum(sum( var1tmp ,3),4);

%var2
var2 = zeros(p,pinp);
for ii = 1:p;
  for jj = 1:pinp;
    var2(ii,jj) = Erhoinv(ii)*SigD(jj,jj);
  end
end

%var3
var3 = squeeze( sum( var3tmp ,2) ); % i,j now

