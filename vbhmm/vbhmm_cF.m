function [Fv]= vbhmm_cF(datastrings,alphabet,net);

%Scores test sequences under a learnt model.
%
%[Fv] = vbhmm_cF(datastrings,alphabet,net);
%
% datastrings (N by T_n) - cell array of strings (can vary in length)
% alphabet (1 by L) - string of characters
% net - learnt model from using net=vbhmm(...)
%
% returns:
%
% Fv (N by 1) - lower bound log likelihood scores for each sequence
%
%M J Beal 14/04/02

L = size(alphabet,2); % cardinality of alphabet
N = size(datastrings,1); % number of sequences to train on
K = size(net.Wa,1);

% convert data into {1,...,L}^N using alphabet
for n = 1:N,
  T(n) = length(datastrings{n});
  tmp=zeros(1,T(n));
  for j=1:L,
    tmp(1,find(datastrings{n}==alphabet(j)))=j;
  end;
  data{n}=tmp;
end;

astar = exp(digamma(net.Wa) - repmat(digamma(sum(net.Wa,2)),1,K));
bstar = exp(digamma(net.Wb) - repmat(digamma(sum(net.Wb,2)),1,L));
pistar = exp(digamma(net.Wpi) - digamma(sum(net.Wpi)));

% E Step only
[wa wb wpi tmp Fv] = forwback(astar,bstar,pistar,data);
