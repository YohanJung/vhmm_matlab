% function to evaluate a test sequence(s) on a learnt model (net)
%
% Returns the 'net' with the training hidden sequences filled in.
%
% function [F,lnZpn,net] = Ftest(Yn,inpn,net);
function [F,lnZpn,net] = Ftest(Yn,inpn,net,dispopt);

if nargin<4, dispopt=0; end

% Initialise dimensions of inputs/outputs
n = size(Yn,2);			% Yn is a cell array of size ( 1 x #sequences )
Tn = zeros(1,n);		% Tn is a row vector containing length of each seq
for i=1:n
  Tn(i) = size(Yn{i},2);	% sequence length may vary
end
p = size(Yn{1},1);		% take first obs seq to set the observ. dim
pinp = size(inpn{1},1);		% take first input sequence to set input dim

if size(inpn,2)~=n		% then use the first input sequence as master
  mstseq=inpn{1};
  if size(mstseq,2)<max(Tn)
    if dispopt
      fprintf('Error: I will only replicate the first inp seq, and it is not long enough'); return;
    end
  end
  if dispopt
    fprintf('[ replicating first inp seq ]\n');
  end
  for i=1:n
    inpn{i}=mstseq(:,1:Tn(i));
  end
end

if dispopt
  fprintf('\n');
  fprintf('# of test sequences               : %i\n',n);
  fprintf('Min,Max sequence lengths          : %i,%i\n',min(Tn),max(Tn));
  fprintf('Total # of observations           : %i\n',sum(Tn,2));
end

cFbool = 1;

net.hidden = [];  % clear old hidden state states (these were for the training seqs, now for the test seqs)
% VBE step : q(X) is Gaussian
for i=1:n
  [lnZpTn{i},Xm_a,Xci_a] = forwardpass(Yn{i},inpn{i},net.hyp.X0ci_p,net.hyp.X0m_p,net.exp.A,net.exp.AA,net.exp.AB,net.exp.B,net.exp.BB,net.exp.rho,net.exp.lnrho,net.exp.CrhoC,net.exp.rhoC,net.exp.CrhoD,net.exp.rhoD,net.exp.DrhoD,cFbool);
  if cFbool==1, lnZpn(i) = sum(lnZpTn{i},2); end
  [Xm_b,Xci_b,X0m_b,X0ci_b] = backwardpass(Yn{i},inpn{i},net.exp.A,net.exp.AA,net.exp.B,net.exp.AB,net.exp.CrhoC,net.exp.rhoC,net.exp.CrhoD);
  [net.hidden.Xmn{i},net.hidden.Xcin{i},net.hidden.Upsn{i},net.hidden.X0mn{i},net.hidden.X0cin{i},net.hidden.Ups0n{i}] = getmarginals(Xm_a,Xci_a,Xm_b,Xci_b,net.hyp.X0m_p,net.hyp.X0ci_p,X0m_b,X0ci_b,net.exp.A,net.exp.AA,net.exp.CrhoC);  
end  

[F] = Fcalc(net.hyp.pa,net.hyp.pb,net.hyp.alpha,net.hyp.beta,net.hyp.gamma,net.hyp.delta,net.param.SigA,net.param.SigB,net.param.SigC,net.param.SigD,net.param.Bb,net.param.Db,net.param.SA,net.param.SC,net.param.GA,net.param.MD,net.param.G,lnZpn,Tn);

