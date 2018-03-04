%Function to perform a variational forward filtering pass on one sequence.
%
%v3.0 Matthew J. Beal GCNU 01/03/03

function [lnZp,mua,nua] = forwardpass(y,u,nup,mup,A,AA,AB,B,BB,rho,lnrho,CrhoC,rhoC,CrhoD,rhoD,DrhoD,cFbool);

k = size(mup,1); [p T] = size(y);

lnZp = zeros(1,T);
nua = zeros(k,k,T);
mua = zeros(k,T);

t=1;
  tmp = inv( nup + AA );
  nua(:,:,t) = eye(k,k) + CrhoC - A*tmp*A';
  mua(:,t) = inv(nua(:,:,t))*( (B-CrhoD-A*tmp*AB)*u(:,t) + rhoC'*y(:,t) +A*tmp*nup*mup );
  if cFbool==1
    lnZp(1,t) = -lndet(nup)-lndet(tmp)+lndet(nua(:,:,t))-sum(lnrho,1)+p*log(2*pi) ...
	+u(:,t)'*BB*u(:,t)+mup'*nup*mup ...
	-(nup*mup-AB*u(:,t))'*tmp*(nup*mup-AB*u(:,t)) ...
	+y(:,t)'*rho*y(:,t)-2*y(:,t)'*rhoD*u(:,t)+u(:,t)'*DrhoD*u(:,t)-mua(:,t)'*nua(:,:,t)*mua(:,t);
  end
for t=2:T
  tmp = inv( nua(:,:,t-1) + AA );
  nua(:,:,t) = eye(k,k) + CrhoC - A*tmp*A';
  mua(:,t) = inv(nua(:,:,t))*( (B-CrhoD-A*tmp*AB)*u(:,t) + rhoC'*y(:,t) +A*tmp*nua(:,:,t-1)*mua(:,t-1) );
  if cFbool==1
    lnZp(1,t) = -lndet(nua(:,:,t-1))-lndet(tmp)+lndet(nua(:,:,t))-sum(lnrho,1)+p*log(2*pi) ...
	+u(:,t)'*BB*u(:,t)+mua(:,t-1)'*nua(:,:,t-1)*mua(:,t-1) ...
	-(nua(:,:,t-1)*mua(:,t-1)-AB*u(:,t))'*tmp*(nua(:,:,t-1)*mua(:,t-1)-AB*u(:,t)) ...
	+y(:,t)'*rho*y(:,t)-2*y(:,t)'*rhoD*u(:,t)+u(:,t)'*DrhoD*u(:,t)-mua(:,t)'*nua(:,:,t)*mua(:,t);
  end
end

if cFbool==1
  lnZp = -.5*lnZp;
else
  lnZp = [];
end
