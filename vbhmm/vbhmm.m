function [net]= vbhmm(datastrings,alphabet,K,its,tol,net);

%Variational Bayesian Hidden Markov Model
%
%net = vbhmm(datastrings,alphabet,K,its,tol,net);
%
% datastrings (N by T_n) - cell array of strings (can vary in length)
% alphabet (1 by L) - string of characters
% K - number of states 
% its - maximum number of iterations of VB EM (100)
% tol - termination tolerance, prop. change in per-datum likelihood (0.0001)
% net - previously learnt model (otherwised initialised from prior)
%
%net is a structure consisting of:
%
% net.Wa - state transition Dirichlet counts
% net.Wb - observation emission Dirichlet counts
% net.Wpi - initial state prior Dirichlet counts
% net.F - F learning curve
%
% Iterates until a proportional change < tol in the per-sequence log
% likelihood, or its iterations of VB.
%
%M J Beal 14/04/02

if nargin<6,
  initfromprior = 1;
else
  initfromprior = 0;
end;
if nargin<5,
  tol = 0.0001;
end;
if nargin<4,
  its = 100;
end;

L = size(alphabet,2); % cardinality of alphabet
N = size(datastrings,1); % number of sequences to train on

% convert data into {1,...,L}^N using alphabet
for n = 1:N,
  T(n) = size(datastrings{n},2);
  tmp = zeros(1,T(n));
  for j=1:L,
    tmp(1,find(datastrings{n}==alphabet(j)))=j;
  end;
  data{n} = tmp;
end;

fprintf('\n********************************************************************\n');
fprintf('Training %i sequences of maximum length %i from an alphabet of size %i\n',N,max(T),L);
fprintf('Variational Bayesian HMM with %i hidden states\n',K);
fprintf('********************************************************************\n\n');
total_length = sum(T);

% Initialise the hyperparameters
alphaa=1;
alphab=1;
alphapi=1;
% Initialise the pseudo-counts
ua = ones(1,K)*(alphaa/K);
ub = ones(1,L)*(alphab/L);
upi = ones(1,K)*(alphapi/K);
% Pick an HMM from the prior to initialize the counts
wa = []; wb = [];
for k=1:K % loop over hidden states
  wa(k,:) = dirrnd(ua,1)*total_length;
  wb(k,:) = dirrnd(ub,1)*total_length;
end
wpi = dirrnd(upi,1)*N;

Fold = -Inf; ntol = tol*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=1:its

  if (initfromprior==0 & it==1)
    Wa=hmm.Wa;
    Wb=hmm.Wb;
    Wpi=hmm.Wpi;
    disp('initing');
  else
    % M Step
    Wa = wa + repmat(ua,[K 1]); % posterior is data counts plus prior.
    Wb = wb + repmat(ub,[K 1]);
    Wpi = wpi + upi;
  end
  
  astar = exp(  digamma(Wa) - repmat( digamma(sum(Wa,2)) ,[1 K])  );
  bstar = exp(  digamma(Wb) - repmat( digamma(sum(Wb,2)) ,[1 L])  );
  pistar = exp(  digamma(Wpi) - digamma(sum(Wpi,2))  );

  % E Step
  [wa wb wpi lnZ(it) lnZv] = forwback(astar,bstar,pistar,data);
  
  % Compute F, straight after E Step.
  Fa(it)=0; Fb(it)=0; Fpi(it)=0;
  for kk = 1:K
    Fa(it) = Fa(it) - kldirichlet(Wa(kk,:),ua);
    Fb(it) = Fb(it) - kldirichlet(Wb(kk,:),ub);
  end
  
  Fpi(it) = - kldirichlet(Wpi,upi);
  
  F(it) = Fa(it)+Fb(it)+Fpi(it)+lnZ(it);

  if it == 1
    fprintf('It:%3i \tFa:%3.3f \tFb:%3.3f \tFpi:%3.3f \tFy:%3.3f \tF:%3.3f\n',it,Fa(it),Fb(it),Fpi(it),lnZ(it),F(it));
  else
    fprintf('It:%3i \tFa:%3.3f \tFb:%3.3f \tFpi:%3.3f \tFy:%3.3f \tF:%3.3f \tdF:%3.3f\n',it,Fa(it),Fb(it),Fpi(it),lnZ(it),F(it),F(it)-Fold);
  end
  
  Fold = F(it);
  if (it>2) 
    if (F(it)<(F(it-1) - 1e-6))    
        fprintf('violation');
    elseif ((F(it)-F(2))<(1 + ntol)*(F(it-1)-F(2)) || ~isfinite(F(it))) 
      fprintf('\nconverged\nend\n');    break;
    end;
  end;
end;

net.Wa = Wa;
net.Wb = Wb;
net.Wpi = Wpi;
net.F = [[1:size(F,2)]' Fa' Fb' Fpi' lnZ' F'];
