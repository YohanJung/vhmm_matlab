%Function to perform a variational backward smoothing pass on one sequence.
%
%v3.0 Matthew J. Beal GCNU 01/03/03

function [Xm_b,Xci_b,X0m_b,X0ci_b] = backwardpass(Y,inp,A,AA,B,AB,CrhoC,rhoC,CrhoD);

k = size(A,1);
[p T] = size(Y);

Xci_b = zeros(k,k,T);
Xm_b  = zeros(k,T);
nus = zeros(k,k);

% beta_T is set to 1, and so it not represented here.  However we leave
% Xci_b(:,:,T) as zero (as initialised) which emulates nicely the same effect
% (in the marginal calculations) as the setting of beta_T=1;

t = T;
  nus            = eye(k,k)+CrhoC;
  Xci_b(:,:,t-1) = AA-A'*inv(nus)*A;
  Xm_b(:,t-1)    = inv(Xci_b(:,:,t-1))*( -AB*inp(:,t) + A'*inv(nus)*( B*inp(:,t) + rhoC'*Y(:,t) - CrhoD*inp(:,t) ) );
for t = [(T-1):-1:2]
  nus            = eye(k,k)+CrhoC+Xci_b(:,:,t);
  Xci_b(:,:,t-1) = AA-A'*inv(nus)*A;
  Xm_b(:,t-1)    = inv(Xci_b(:,:,t-1))*( -AB*inp(:,t) + A'*inv(nus)*( B*inp(:,t) + rhoC'*Y(:,t) - CrhoD*inp(:,t) + Xci_b(:,:,t)*Xm_b(:,t) ) );
end
t = 1;
  nus            = eye(k,k)+CrhoC+Xci_b(:,:,t);
  X0ci_b         = AA-A'*inv(nus)*A;
  X0m_b          = inv(X0ci_b          )*( -AB*inp(:,t) + A'*inv(nus)*( B*inp(:,t) + rhoC'*Y(:,t) - CrhoD*inp(:,t) + Xci_b(:,:,t)*Xm_b(:,t) ) );
