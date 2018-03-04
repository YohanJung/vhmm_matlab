function  [wa, wb, wpi, lnZ, lnZv] = forwback(A,B,Pi,data);

%Performs forward-backward message passing for HMMs
%
%[wa, wb, wpi, lnZ, lnZv] = forwback(A,B,Pi,data);
%
% A (K by K) - state transition probabilities
% B (P by K) - observation emission probabilities
% Pi (1 by K) - initial state prior
% Note these probabilities can be sub-normalised.
%
%The E step returns:
%
% wa (K by K) - state transition counts (from all data)
% wb (P by K) - observation emission counts (from all data)
% wpi (1 by K) - initial state prior counts (from all data)
% lnZ (1 by 1) - log likelihood of all data
% lnZv (N by 1) - log likelihood of each data string
%
%M J Beal 13/04/02

[K L] = size(B);
N = size(data,2);
for n = 1:N,
  T(n) = length(data{n});
end;

Xi = zeros(K,K);
Gammak = zeros(K,L);
GammaInit = zeros(1,K);
lnZv = zeros(N,1);

for n=1:N,

  Gamma = zeros(T(n),K);
  alpha = zeros(T(n),K);
  beta = zeros(T(n),K);
  scale = zeros(1,T(n));
  
  % Pick out the likelihood of each symbol in the sequence
  Obs = B(:,data{n})'; 
  
  % Forward pass (with scaling)
  alpha(1,:) = Pi.*Obs(1,:);
  scale(1) = sum(alpha(1,:));
  alpha(1,:) = alpha(1,:)/scale(1);
  for t=2:T(n)
    alpha(t,:) = (alpha(t-1,:)*A).*Obs(t,:);
    scale(t) = sum(alpha(t,:));
    alpha(t,:) = alpha(t,:)/scale(t);
  end;
  
  % Backward pass (with scaling)
  beta(T(n),:) = ones(1,K)/scale(T(n));
  for t=T(n)-1:-1:1
    beta(t,:) = (beta(t+1,:).*Obs(t+1,:))*A'/scale(t); 
  end;
  
  % Another pass gives us the joint probabilities
  for t=1:T(n)-1
    Xi=Xi+A.*(alpha(t,:)'*(beta(t+1,:).*Obs(t+1,:)));
  end
  
  % Compute Gamma
  Gamma = alpha.*beta;
  Gamma = Gamma./repmat(sum(Gamma,2),1,K);
  
  % Compute the sums of Gamma conditioned on k
  for t = 1:T(n)
    Gammak(:,data{n}(t)) = Gammak(:,data{n}(t)) + Gamma(t,:)';
  end;
  
  GammaInit = GammaInit + Gamma(1,:); 

  lnZv(n) = sum(log(scale));
  
end; % for n

lnZ = sum(lnZv,1);
wa = Xi;
wb = Gammak;
wpi = GammaInit;
