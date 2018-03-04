clear all; close all;

set(0,'DefaultAxesPosition',[.1 .1 .8 .8]);
set(0, 'DefaultFigureMenuBar', 'none');
colordef none;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINUS  Demonstration script for vblds model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True model has k = 2 dimensional hidden state-space, and
% let's generate a couple of p = 4 dimensional sequences,
% driven by an input sequence of dimension pinp = 3;
%
% Note, even though the input sequence has 3 dimensions, only the 
% first two will be used to modulate the output (i.e. last column 
% of D contains zeros).

k = 2; p = 4; pinp = 3; % T = 200;
T = {100,300}; N = size(T,2);

% Now, generate the parameters of the true model randomly.

% Let's generate a random A matrix k by k
A = randn(2,2)/2;
% C - p by k output (embedding) matrix,
% drawn from bimodal -2,+2 Gaussians, s.d. 1/2.
sd = sqrt(1/4); C = ( (randn(p,k)+2/sd).*(-1).^(randn(p,k)>0) )*sd;
% Q - wlog is the identity
Q = eye(k,k);
% R - identity, i.e. about 1/5 sd noise
R = eye(p,p);
% Initial (prior) variance and state as usual
V0=Q; x0=zeros(k,1);
% B - k by u hidden state input loadings,
% drawn from indep +/- U(-1,1);
%B = 10*(rand(k,u)-.5);
B = zeros(k,pinp);	% let's say the inputs do not shift the hidden dynamics.
% D - p by u output space input loadings,
% drawn from indep +/- U(-1,1);
D = [20*(rand(p,pinp-1)-.5) zeros(p,1)];

for n = 1:N
  Tn = T{n};
  twoperiods = 100;
  % inp - u by Tn input sequence
  ti = [1:Tn]/twoperiods*4*pi;
  inpn{n} = [cos(ti); sin(ti); rand(1,Tn)*2-1]; % assumes u = 2;

  [Yn{n},xtrue{n}] = ldsinp(Tn,A,C,Q,R,x0,V0,B,D,inpn{n});
end

figure(1);
for n = 1:N
  subplot(3,N,N*0+n); plot(inpn{n}','o-','markersize',4); ylabel('Inputs u_{1:T}'); set(gca,'xticklabel',[]); title(['sequence ' num2str(n)]);
  subplot(3,N,N*1+n); plot(xtrue{n}','o-','markersize',4); if n==1, ylabel('True x_{1:T}'); end; set(gca,'xticklabel',[]);
  subplot(3,N,N*2+n); plot(Yn{n}','o-','markersize',4); if n==1, ylabel('Observations y_{1:T}'); end; xlabel('Time t');
end
  disp(['Figure 1 shows generation process for ' num2str(N) ' sequences:']);
  disp(['  a ' num2str(pinp) '-dim input u_{1:T} - subplot 1,']);
  disp(['  driving a ' num2str(k) '-dim hidden state dynamics x_{1:T} - subplot 2,']);
  disp(['  producing a ' num2str(p) '-dim output y_{1:T} - subplot 3.']);
  
disp(' '); disp('<pause>'); pause

figure(2);
  subplot(221); hinton(A,['true A  ' num2str(max(max(A)),'%.4f')],'standard');
  subplot(222); hinton(B+eps,['true B  ' num2str(max(max(B)),'%.4f')],'standard');
  subplot(223); hinton(C,['true C  ' num2str(max(max(C)),'%.4f')],'standard');
  subplot(224); hinton(D,['true D  ' num2str(max(max(D)),'%.4f')],'standard');
  
  disp(' ');
  disp('Figure 2 shows the true A, B, C and D matrices used in Hinton diagram form:');
  disp('  Each small square represents an entry, with white positive and black negative,');
  disp('  whose area is proportional to its (log) magnitude.');
  disp('  The number above each diagram is the largest value of the weights in each matrix.');

disp(' '); disp('<pause>'); pause

kmodel = 4;

% Prior mean of matrice A,B,C and D
muA= zeros(kmodel,kmodel); muB = zeros(kmodel,pinp); muC= zeros(p,kmodel); muD = zeros(p,pinp);
PriMu{1} = muA;   PriMu{2} = muB;  PriMu{3} = muC; PriMu{4} = muD;

disp(' ');
disp('We''ll now run a Variational Bayesian State-Space Model, with');
disp(['with a state space of dimension ' num2str(kmodel) '.']);

net = vssminpn(Yn,inpn,kmodel,200,3,1,PriMu,[0 0 0]); % Note last arguments 'net' not specified.


%1./[net.hyp.alpha net.hyp.gamma ; NaN NaN ; net.hyp.beta net.hyp.delta]

figure(5);
  subplot(221); hinton(net.exp.A,['post. mean of A  ' num2str(max(max(net.exp.A)),'%.4f')],'standard');
  subplot(222); hinton(net.exp.B+eps,['post. mean of B  ' num2str(max(max(net.exp.B)),'%.4f')],'standard');
  subplot(223); hinton(net.exp.C,['post. mean of C  ' num2str(max(max(net.exp.C)),'%.4f')],'standard');
  subplot(224); hinton(net.exp.D,['post. mean of D  ' num2str(max(max(net.exp.D)),'%.4f')],'standard');
  
figure(6);
for n = 1:N
  subplot(2,N,N*0+n)
    plot(xtrue{n}','o-','markersize',4); ylabel('True x_{1:T}'); set(gca,'xticklabel',[]); title(['sequence ' num2str(n)]);
  subplot(2,N,N*1+n)
    plot(net.hidden.Xmn{n}','o-','markersize',4); ylabel('Inferred x_{1:T}'); set(gca,'xticklabel',[]); xlabel('time');
end

disp(' ');
disp('To run more iterations, simply supply the network ''net'' that we just learnt.');
disp(' '); disp('<pause>'); pause

net = vssminpn(Yn,inpn,kmodel,200,3,1,PriMu,[0 0 0],net); % Note last two arguments left, they will default.

figure(9);
  subplot(221); hinton(net.exp.A,['post. mean of A  ' num2str(max(max(net.exp.A)),'%.4f')],'standard');
  subplot(222); hinton(net.exp.B+eps,['post. mean of B  ' num2str(max(max(net.exp.B)),'%.4f')],'standard');
  subplot(223); hinton(net.exp.C,['post. mean of C  ' num2str(max(max(net.exp.C)),'%.4f')],'standard');
  subplot(224); hinton(net.exp.D,['post. mean of D  ' num2str(max(max(net.exp.D)),'%.4f')],'standard');
  
figure(10);
for n = 1:N
  subplot(2,N,N*0+n)
    plot(xtrue{n}','o-','markersize',4); ylabel('True x_{1:T}'); set(gca,'xticklabel',[]); title(['sequence ' num2str(n)]);
  subplot(2,N,N*1+n)
    plot(net.hidden.Xmn{n}','o-','markersize',4); ylabel('Inferred x_{1:T}'); set(gca,'xticklabel',[]); xlabel('time');
end
