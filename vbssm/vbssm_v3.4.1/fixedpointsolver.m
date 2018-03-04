%caution required in this routine, since digamma() and trigamma() are not defined
%for arguments<0.  Therefore, should a Newton-Raphson step take us out of their
%domain, then we reset a to mina, and carry on the optimisation.
%
%This could be done more accurately using exponential domain (see appendix in thesis).
%
%v3.0 Matthew J. Beal GCNU 01/03/03

function [a,b] = fixedpointsolver(a,b,c,d);

e=log(d)-c;

da=Inf; mina = .005;

maxit=100; it=0;
while abs(da)>1e-6
  it=it+1;
  da = ( digamma(a)-log(a)+e )/( a*trigamma(a)-1 );
  %fprintf('da: %3.4f\n',da);
  %a = a*(1-da);
  a = a*exp(-da);
  if a<mina
    a = max(a,mina);
    fprintf('-');
    da=Inf;
  end
  if it>maxit
    fprintf('*'); break
  end
end

b = a/d;
