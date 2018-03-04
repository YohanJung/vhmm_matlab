%Function to calculate single node marginals and cross-time (pairwise) marginals.
%
%v3.0 Matthew J. Beal GCNU 01/03/03

function [Xm_Y,Xci_Y,Ups_Y,X0m_Y,X0ci_Y,Ups0_Y] = getmarginals(Xm_a,Xci_a,Xm_b,Xci_b,X0m_p,X0ci_p,X0m_b,X0ci_b,A,AA,CrhoC)

[k T] = size(Xm_a);
Xci_Y = zeros(k,k,T);
Xm_Y = zeros(k,T);
Ups_Y = zeros(k,k,T-1);

  t = 0;
  X0ci_Y       =                           X0ci_p           +       X0ci_b            ;
  X0m_Y        =       inv(X0ci_Y)*(       X0ci_p*X0m_p     +       X0ci_b*X0m_b     );
for t = 1:T
  Xci_Y(:,:,t) =                     Xci_a(:,:,t)           + Xci_b(:,:,t)            ;
  Xm_Y(:,t)    = inv(Xci_Y(:,:,t))*( Xci_a(:,:,t)*Xm_a(:,t) + Xci_b(:,:,t)*Xm_b(:,t) );
end


  Ups0_Y       = inv(X0ci_p      +AA) *A' *inv( eye(k,k) + CrhoC + Xci_b(:,:,1)   - A*inv(X0ci_p      +AA)*A' );
for t = 1:(T-1)
  Ups_Y(:,:,t) = inv(Xci_a(:,:,t)+AA) *A' *inv( eye(k,k) + CrhoC + Xci_b(:,:,t+1) - A*inv(Xci_a(:,:,t)+AA)*A' );
end
