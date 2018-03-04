% function [lik, likv]=dhmm_cl(X,alphabet,E,P,Pi);    ...also...
%
% function [lik, likv]=maphmm_cl(X,alphabet,E,P,Pi);
% 
% simple Hidden Markov Model - variable lengths and discrete symbols
% calculate likelihoods for X
%
% X - cell array of strings
% alphabet - string of characters
% E - observation emission probabilities
% P - state transition probabilities
% Pi - initial state prior probabilities
%
% lik - total log likelihood 
% likv - vector of log likelihoods for each sequence

function [lik, likv]=dhmm_cl(X,alphabet,E,P,Pi);

LA = length(alphabet); 
epsi=1e-100;

% number of sequences
N=length(X);

% length of each sequence
T=ones(1,N);
for n=1:N,
  T(n)=length(X{n});
end;

TMAX = max(T);

[dummy,K]=size(E);

if (dummy ~= LA | size(P) ~=[K K] | size(Pi) ~= [1 K]) 
   error('dhmm_cl: error: incorrect sizes');
end;		   

B=zeros(TMAX,K);

LL=[];
likv=zeros(N,1);
lik=0;

%%%% FORWARD 

for n=1:N
    alpha=zeros(T(n),K);

    % Inital values of B = Prob(output|s_i), given data X

    Xcurrent=X{n};
    
    for i=1:T(n)
      m = findstr(alphabet,Xcurrent(i));
      if (m == 0)
	fprintf('symbol not found\n');  
	return;
      end
      B(i,:) = E(m,:);
    end;

    scale=zeros(T(n),1);
    alpha(1,:)=Pi(:)'.*B(1,:)+epsi;
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for i=2:T(n)
      alpha(i,:)=(alpha(i-1,:)*P).*B(i,:)+epsi;
      scale(i)=sum(alpha(i,:));
      alpha(i,:)=alpha(i,:)/scale(i);
    end;

    likv(n)=sum(log(scale(1:T(n),:)));

end;

lik=sum(likv);




