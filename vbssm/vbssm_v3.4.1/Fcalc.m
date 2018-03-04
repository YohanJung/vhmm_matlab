% Function to calculate the lower bound F, for the variational Bayesian state-space model.
%
% Broken down into 6 parts.  These should be obvious with reference to thesis.
%
% Matthew J. Beal GCNU 04/01/06

function [F] = Fcalc(pa,pb,alpha,beta,gamma,delta,SigA,SigB,SigC,SigD,Bb,Db,SA,SC,GA,MD,G,lnZpn,Tn,PriMu) 

muA = PriMu{1};
muB = PriMu{2};
muC = PriMu{3};
muD = PriMu{4};

F = zeros(1,6);
k = size(alpha,1);
pinp = size(beta,1);
p = size(G,1);

T = sum(Tn,2); 

F(1) = -( -k/2*( sum(log(beta))+lndet(SigB) ) ...
	  + k/2*trace( diag(beta)*SigB-eye(pinp,pinp) ) ...
      + 1/2*trace( diag(beta)*(SigB*Bb-muB')*(SigB*Bb-muB')') ); 
      
F(2) = -( -k/2*( sum(log(alpha))+lndet(SigA) )...
	  + k/2*trace( diag(alpha)*SigA-eye(k,k) ) ...
      + k/2*trace( diag(alpha)*muA'*muA ) ...                      
	  + k/2*trace( diag(alpha)*SigA*GA*SigB*GA'*SigA ) ...
      +1/2*trace( diag(alpha)*SigA*(SA+diag(alpha)*muA'-GA*SigB*Bb)*(SA+diag(alpha)*muA'-GA*SigB*Bb)'*SigA ) ...
      - trace(diag(alpha)*SigA*(SA+diag(alpha)*muA'-GA*SigB*Bb)*muA ) );

 
    
F(3) = - klgamma(pa+sum(Tn,2)/2,pb+G'/2,pa,pb);

F(4) = -( -p/2*( sum(log(delta))+lndet(SigD) ) ...
	  + p/2*trace( diag(delta)*SigD-eye(pinp,pinp) ) ...
      + 1/2*trace( diag(delta)*(SigD*Db-muD')*(pa+sum(Tn,2)/2)*diag(1./(pb+G/2))*(SigD*Db-muD')' ) );
   
  
F(5) = -( -p/2*( sum(log(gamma))+lndet(SigC) ) ...
	  + p/2*trace( diag(gamma)*SigC-eye(k,k) ) ...
      + p/2*trace( diag(gamma)*muC'*(pa+sum(Tn,2)/2)*diag(1./(pb+G/2))*muC ) ... 
	  + p/2*trace( diag(gamma)*SigC*MD*SigD*MD'*SigC ) ...
	  + 1/2*trace( diag(gamma)*SigC*(SC+diag(gamma)*muC'-MD*SigD*Db)*(pa+sum(Tn,2)/2)*diag(1./(pb+G/2))*(SC+diag(gamma)*muC'-MD*SigD*Db)'*SigC ) ...
      - trace( diag(gamma)*SigC*(SC+diag(gamma)*muC'-MD*SigD*Db)*(pa+sum(Tn,2)/2)*diag(1./(pb+G/2))*muC )  );

  
F(6) = sum(lnZpn,2);
