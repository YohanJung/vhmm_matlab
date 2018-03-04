%Generates data from a state-space model (LDS) with inputs if specified:
%
%[y,x] = lds(T,A,C,Q,R,x0,V0,B,D,inp)
%
%  x_t = A x_{t-1} + B u_t + N(0,Q)
%  y_t = C x_t + D u_t + N(0,R)
%
% x,y - hidden and observation sequences.
%
% T - is length of sequence (if vector multiple sequences generated).
% A (k by k) - hidden state dynamics.
% C (p by k) - output (embedding) matrix.
% Q (k by k) - hidden state noise covariance.
% R (p by p) - output noise covariance.
% x0 (k by 1) - initial state mean (zero).
% V0 (k by k) - initial state covariance (Q).
%
%If using inputs:
%
% B (k by u) - hidden state input loadings.
% D (p by u) - output space input loadings.
% inp (u by T_u) - input sequence; T_u must be > max(T);
%
%If T a vector, x and y become cell arrays.
%
%Matthew J. Beal GCNU 04/04/02
%(Input-version of Sam Roweis' lds.m code.)

function [y,x] = ldsinp(T,A,C,Q,R,x0,V0,B,D,inp)

if(nargin<7) V0=[]; end
if(nargin<6) x0=[]; end

[p k] = size(C); [k1 k2] = size(A); assert(all([k1 k2]==k));

x0=x0(:); % makes into a col vector

% tricks to speed up input
if(all(size(Q)==[1 1])) Q=Q*eye(k,k); end
if(all(size(R)==[1 1])) R=R*eye(p,p); end
if(all(size(V0)==[1 1])) V0=V0*eye(k,k); end
if(all(size(Q)==[1 k]) | all(size(Q)==[k,1]) ) Q=diag(Q); end
if(all(size(R)==[1 p]) | all(size(R)==[p,1]) ) R=diag(R); end
if(all(size(V0)==[1 k]) | all(size(V0)==[k,1]) ) V0=diag(V0); end

if(isempty(Q)) Q=eye(k,k); end
if(isempty(R)) R=eye(p,p); end
if(isempty(V0)) V0=Q; end
if(isempty(x0)) x0=zeros(k,1); end

if(nargin<10) inp=zeros(1,max(T)); end
if(nargin<9) D=zeros(p,1); end
if(nargin<8) B=zeros(k,1); end
[u Tu] = size(inp); [k3 u2] = size(B); [p1 u1] = size(D); assert(all([k3==k p1==p u1==u u2==u Tu>=max(T)]));

assert(all(size(Q) ==[k k])); assert(all(size(R) ==[p p]));
assert(all(size(x0)==[k 1]));  assert(all(size(V0)==[k k]));


numseqs=length(T);
if(numseqs==1) 
  [y,x]=genseq(T,A,C,Q,R,x0,V0,B,D,inp);
else
  y=cell(numseqs,1); x=cell(numseqs,1);
  for numseq=1:numseqs
    [Y{numseq},X{numseq}]=genseq(T(numseq),A,C,Q,R,x0,V0,B,D,inp);
  end
end


function [y,x] = genseq(T,A,C,Q,R,x0,V0,B,D,inp)

% initialize
[p k] = size(C);
y = zeros(p,T);
x = zeros(k,T);

sqQ = sqrtm(Q); 

% by resetting the random seed in this way, we can be sure of generating
% incremental datasets. M J Beal

seed=rand*1000;
resetr(seed);
% generate hidden states
x(:,1)=sqrtm(V0)*randn(k,1)+x0+B*inp(:,1);

for t=2:T
  %x(:,t) = A*x(:,t-1)+B*inp(:,1)+sqQ*randn(k,1);
  x(:,t) = A*x(:,t-1)+B*inp(:,t)+sqQ*randn(k,1);
end

resetr(seed+seed);
% generate observations
y = C*x + D*inp(:,1:T) + sqrtm(R)*randn(p,T);

% the rand operation goes down rows first then across time steps, so this is
% fine too for incremental datasets.