%|              Variational Bayesian State-Space models                 |
%|                                                                      |
%|      Copyright Matthew J. Beal 22/05/01 version 3.0 (11/08/03)       |
%|              http://www.cs.toronto.edu/~beal/software                |
%
%USAGE:  net = vssminpn(Yn,inpn,k,its,dispopt,cFbool,PriMu,hypit,net);
%
%-Reqrd-
% Yn      N obs sequences    cell (1 by N), each cell is matrix (p by Tn)
% inpn    input sequence(s)  as for Yn, but can also be just one sequence
%
%-Optnl-  description (df default setting)              (possible values)
% k       maximum hidden state dimension (df p)              (at least 1)
% its     maximum it#,iterations (df 100)                    (at least 1)
% dispopt verbosity of output (0=nothing, df 1)                 (0,1,2,3)
% cFbool  calc F @each it (df 0, but F always calc'd at last it)    (0,1)
% PriMu   prior mean matrix of A,B,C,D
% hypit   it# to optimise dyn,oput,&other hyps (df [10 20 30])   (1 by 3)
% net     a previous output (overwrites k setting)     (matlab structure)
%
%Returns 'net', a structure consisting of:
%
% Optim. hyperparams: net.hyp.{pa,pb,alpha,beta,gamma,delta,X0m_p,X0ci_p}
% Sufficient statistics of the parameter posterior:  net.param.{SigA,...}
% Expected natural parameters:                    net.exp.{A,B,AA,AB,...}
% Sufficient stastistics of the hidden state:
%                             net.hidden.{X0mn,Xmn,X0cin,Xcin,Ups0n,Upsn}
% Histories/trajectories of the lower bound and hyperparameters:
%                                            net.hist.{F,pa,pb,alpha,...}
% E.g. net.exp.B is the mean dynamics matrix, and
%      net.param.SigB is the covariance of each of its rows.
% Note: not the same for matrix A! (see ch.5 of my thesis)
%
%v3.4.1 Matthew J. Beal GCNU 04/01/06

function net = vssminpn(Yn,inpn,k,its,dispopt,cFbool,PriMu,hypit,net);

% Initialise dimensions of inputs/outputs
n = size(Yn,2);			% Yn is a cell array of size ( 1 x #sequences )
Tn = zeros(1,n);		% Tn is a row vector containing length of each seq
for i=1:n
  Tn(i) = size(Yn{i},2);	% sequence length may vary
end
p = size(Yn{1},1);		% take first obs seq to set the observ. dim
pinp = size(inpn{1},1);		% take first input sequence to set input dim

% If just one input sequence is specified, replicate it to be used as the
% input sequence for every sequence.
if size(inpn,2)~=n		% then use the first input sequence as master
  mstseq=inpn{1};
  if size(mstseq,2)<max(Tn)
    fprintf('Error: I will only replicate the first inp seq, and it is not long enough'); return;
  end
  fprintf('[ replicating first inp seq ]\n');
  for i=1:n
    inpn{i}=mstseq(:,1:Tn(i));
  end
end

% Compute some statistics of the data & inputs (no need to put in VBEM loop)
Ydd = zeros(p,p); Udd=zeros(pinp,pinp); UY=zeros(pinp,p);
for i=1:n
  Ydd = Ydd+Yn{i}*Yn{i}';
  Udd = Udd+inpn{i}*inpn{i}';
  UY = UY+inpn{i}*Yn{i}';
end

if nargin<3, k = p; end		% if k not specified, set equal to dimensions of output
				%   this irrelevant if 'net' is provided in function call
if nargin<4, its = 100; end	% if no limit specified then set to 100 iterations
if nargin<5, dispopt=2; end	% degree of verbosity set to just print lower bound progress (no plots)
if nargin<6, cFbool=0; end      % default is to leave F uncalculated (faster)
if nargin<7,               % if no prior mean for matrice A,B,C,D specified then set to zero mean. 
    muA=zeros(k,k);
    muB=zeros(k,pinp);
    muC=zeros(p,k);
    muD=zeros(p,pinp);
else 
    muA = PriMu{1};
    muB = PriMu{2};
    muC = PriMu{3};
    muD = PriMu{4};
end

if nargin<8, 
  dynhypit=10;			% it# after which dynamics hyps optimised
  outhypit=20;			% it# after which output hyps optimised
  etchypit=30;			% it# after which other (input & noise hyps) are optimised
else
  dynhypit=hypit(1); outhypit=hypit(2); etchypit=hypit(3); % user-specified schedule
end

if nargin<9, % Set up all hyperparameters and initialisation right now
  if dispopt>=1; fprintf('Initialising hyps,'); end
  X0m_p = zeros(k,1);		% time zero mean of X set to zero
  X0ci_p = eye(k,k);		% time zero precision X set to identity
  pa = 1; pb = 1;		% shape and precision of noise precision prior
  alpha = 1*ones(k,1);		% precision vector, each element for a COL of A
  beta = 1*ones(pinp,1);	% precision vector, each element for a COL of B
  gamma = 1*ones(k,1);		% pseudo-precision vector, each element for a COL of C
  delta = 1*ones(pinp,1);	% pseudo-precision vector, each element for a COL of D
  
  
  % Initialise a randomised hidden state
  %   This is somehow preferable to initialising a set of expected natural parameters
  %   from their prior, which technically speaking is the more correct thing to do --- Sorry!
  if dispopt>=1; fprintf(' and randomly instantiating hidden state sequences'); end
  for i=1:n
    X0mn{i} = X0m_p;				% hidden state means
    Xmn{i} = randn(k,Tn(i));
    X0cin{i} = X0ci_p;				% hidden state inverse covariances
    Xcin{i} = repmat(eye(k,k),[1 1 Tn(i)]);
    Ups0n{i} = eye(k,k);			% hidden state cross-correlations
    Upsn{i} = repmat(eye(k,k),[1 1 Tn(i)]);
  end
  % Set up hyperparameter histories
  pahist = pa; pbhist = pb; alphahist = alpha; betahist = beta; gammahist = gamma; deltahist = delta; X0m_phist = X0m_p; X0ci_phist = X0ci_p;
  hist = -Inf*ones(1,7);

else				% Use the net provided in the function call

  if dispopt>=1; fprintf('Extracting hyperparameters, their histories, and parameter posterior  from previously learnt network'); end
  % Extract the hyperparameters...
  pa	= net.hyp.pa;
  pb	= net.hyp.pb;
  alpha	= net.hyp.alpha; k=size(net.hyp.alpha,1);
  beta	= net.hyp.beta;
  gamma	= net.hyp.gamma;
  delta	= net.hyp.delta;
  X0m_p	= net.hyp.X0m_p;
  X0ci_p= net.hyp.X0ci_p;
  % Extract also the expected natural parameters...
  A	= net.exp.A;
  B	= net.exp.B;
  AA	= net.exp.AA;
  AB	= net.exp.AB;
  BB	= net.exp.BB;
  rho	= net.exp.rho;
  lnrho	= net.exp.lnrho;
  CrhoC	= net.exp.CrhoC;
  rhoC	= net.exp.rhoC;
  CrhoD	= net.exp.CrhoD;
  rhoD	= net.exp.rhoD;
  DrhoD	= net.exp.DrhoD;
  % Finally now infer the hidden state sequence in the 0th VBEM step.
  %   NOTA BENE: this way it is possible to use the posterior over parameters 
  %   for training set A, and then continue training FROM THAT POINT on a second 
  %   training set B (which may have sequences of different length).
  % VBE step : q(X) is Gaussian
  for i=1:n
    [lnZpTn{i},Xm_a,Xci_a] = forwardpass(Yn{i},inpn{i},X0ci_p,X0m_p,A,AA,AB,B,BB,rho,lnrho,CrhoC,rhoC,CrhoD,rhoD,DrhoD,cFbool);
    [Xm_b,Xci_b,X0m_b,X0ci_b] = backwardpass(Yn{i},inpn{i},A,AA,B,AB,CrhoC,rhoC,CrhoD);
    [Xmn{i},Xcin{i},Upsn{i},X0mn{i},X0cin{i},Ups0n{i}] = getmarginals(Xm_a,Xci_a,Xm_b,Xci_b,X0m_p,X0ci_p,X0m_b,X0ci_b,A,AA,CrhoC);  
  end  
  % Gather up previous hyperparameter histories
  pahist = net.hist.pa;
  pbhist = net.hist.pb;
  alphahist = net.hist.alpha;
  betahist = net.hist.beta;
  gammahist = net.hist.gamma;
  deltahist = net.hist.delta;
  X0m_phist = net.hist.X0m_p;
  X0ci_phist = net.hist.X0ci_p;
  hist = net.hist.F;

end

% Display some summary information on the data sequences and proposed model size
if dispopt>=1;
  fprintf('.\nSummary:');
  fprintf(' # of sequences           : %i\n',n);
  fprintf('         Min,Max sequence lengths : %i,%i\n',min(Tn),max(Tn));
  fprintf('         Total # of observations  : %i\n',sum(Tn,2));
  fprintf('         # hidden states in model : %i',k);
  if ~cFbool, fprintf('\nSkipping F calculations (fast, but no history)'); end
end

% Set up figure plots:
if dispopt>=3,
  hf = figure; clf
  subplot(211);
    hFpart = plot(zeros(2,6)); title('F constituents');
    legend('B','A|B','rho','D|rho','C|rho,D','F(Y)');
  subplot(212);
    hFdiff = plot(zeros(1,1)); title('log rate of change of F');
    legend('log \DeltaF(Y)');
  hf2 = figure; clf; 
  subplot(231); halpha = plot(zeros(2,k)); title('prior log var A_{|.}');
  subplot(232); hbeta = plot(zeros(2,pinp)); title('prior log var B_{|.}');
  subplot(234); hgamma = plot(zeros(2,k)); title('prior log var C_{|.}');
  subplot(235); hdelta = plot(zeros(2,pinp)); title('prior log var D_{|.}');
  subplot(233); hrhomean = plot(zeros(2,p)); title('prior \rho mean');
  subplot(236); hrhovar = plot(zeros(2,p)); title('prior \rho log var');
  set(gcf,'doublebuffer','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN VBEM iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

it = 0; dF=Inf; tic;
while 1 % note there is a loop break after the F-calculation, of the sort: ~( (it<its) & (dF>1e-14) )

  it = it+1; if (it==4) & (dispopt>=1), fprintf('\nApproxTime for %4i iterations     : %2.1f mins (''.''=50its)',its,toc/3*(its-3)/60); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % UPDATE A,B
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % VBM step (dynamics model)
  % Compute parameters to represent Q(A,B)
  WA = zeros(k,k);
  SA = zeros(k,k);
  GA = zeros(k,pinp);
  MB = zeros(k,pinp);
  for i=1:n
      WA = WA+inv(double(X0cin{i}))+ double(X0mn{i})*double(X0mn{i}')+ double(Xmn{i}(:,1:(Tn(i)-1)))*double(Xmn{i}(:,1:(Tn(i)-1))');
      for t = 1:(Tn(i)-1);
          WA = WA+inv(double(Xcin{i}(:,:,t)));
      end
      SA = SA+ Ups0n{i}+double(X0mn{i})*double(Xmn{i}(:,1)') + sum(Upsn{i}(:,:,1:(Tn(i)-1)),3)+double(Xmn{i}(:,1:Tn(i)-1))*double(Xmn{i}(:,2:Tn(i))');

      GA = GA+ double(X0mn{i})*double(inpn{i}(:,1)') + double(Xmn{i}(:,1:Tn(i)-1))*double(inpn{i}(:,2:Tn(i))');
      MB = MB+ double(Xmn{i})*double(inpn{i}');
  end
  
  SigA = inv( WA + diag(alpha) );
  SigB = inv( Udd - GA'*SigA*GA + diag(beta) );
  

  Bb = MB' - GA'*SigA*(SA+diag(alpha)*muA') + diag(beta)*muB'; 

  
  % VBM step compute sufficient statistics (dynamics model)
  B = Bb'*SigB;					% < B >

  Ab = SA - GA*B' + diag(alpha)*muA';     
  
  A = Ab'*SigA;		      % < A > 
  
  AA = A'*A + k*( SigA + SigA*GA*SigB*GA'*SigA );	% < A' A >
  
  AB = SigA*(SA+diag(alpha)*muA')*B - SigA*GA*( B'*B + k*SigB );		% < A' B > 
  BB = B'*B + k*SigB;					% < B' B >
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % UPDATE rho,C,D
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % VBM step (observation model)
  % Compute parameters to represent Q(rho,C,D)
  WC = zeros(k,k);
  SC = zeros(k,p);
  MD = zeros(k,pinp);   % Nota Bene, MD is equivalent to GC in thesis.
  for i=1:n
    WC = WC+Xmn{i}*Xmn{i}';
    for t = 1:Tn(i)
      WC = WC+inv(Xcin{i}(:,:,t));
    end
    SC = SC+ Xmn{i}*Yn{i}';
    MD = MD+ Xmn{i}*inpn{i}';
  end
  SigC = inv( WC + diag(gamma) );
  SigD = inv( Udd - MD'*SigC*MD + diag(delta) );
  
  
  Db = UY - MD'*SigC*(SC+diag(gamma)*muC') + diag(delta)*muD';   

  G = diag( Ydd - (SC+diag(gamma)*muC')'*SigC*(SC+diag(gamma)*muC')...
              -Db'*SigD*Db + muD*diag(delta)*muD' );                  
  
  % VBM step compute sufficient statistics (observation model)
  rho = (pa+sum(Tn,2)/2)*diag(1./(pb+G/2));		% < R^{-1} > matrix
  lnrho = digamma(pa+sum(Tn,2)/2) - log(pb+G/2);	% < -log R_ss > vector
  D = Db'*SigD;						% < D >
  
  Cb = SC - MD*D'+ diag(gamma)*muC';       
  C = Cb'*SigC;				% < C >
  CrhoC = C'*rho*C + p*( SigC + SigC*MD*SigD*MD'*SigC );% < C' R^{-1} C >
  rhoC = rho*C;						% < R^{-1} C >
  CrhoD = SigC*( SC*rho*D - MD*D'*rho*D - p*MD*SigD );	% < C' R^{-1} D >
  rhoD = rho*D;						% < R^{-1} D >
  DrhoD = D'*rho*D + p*SigD;				% < D' R^{-1} D >

  % Is this the last iteration?
  finalit = ~(it<its);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % UPDATE hidden state
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % VBE step : q(X) is Gaussian
  for i=1:n
    [lnZpTn{i},Xm_a,Xci_a] = forwardpass(Yn{i},inpn{i},X0ci_p,X0m_p,A,AA,AB,B,BB,rho,lnrho,CrhoC,rhoC,CrhoD,rhoD,DrhoD,cFbool|finalit);
    if cFbool|finalit , lnZpn(i) = sum(lnZpTn{i},2); end
    [Xm_b,Xci_b,X0m_b,X0ci_b] = backwardpass(Yn{i},inpn{i},A,AA,B,AB,CrhoC,rhoC,CrhoD);
    [Xmn{i},Xcin{i},Upsn{i},X0mn{i},X0cin{i},Ups0n{i}] = getmarginals(Xm_a,Xci_a,Xm_b,Xci_b,X0m_p,X0ci_p,X0m_b,X0ci_b,A,AA,CrhoC);  
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE LOWER BOUND, F
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if cFbool|finalit
    % Calculate F = -kl_A -kl_B -kl_C -kl_D -kl_rho + \sum_{i=1}^n lnZpn_i
    % For this we require parameter sufficient statistics and lnZp
    [F] = Fcalc(pa,pb,alpha,beta,gamma,delta,SigA,SigB,SigC,SigD,Bb,Db,SA,SC,GA,MD,G,lnZpn,Tn,PriMu);
    % book keeping for F constituents
    hist = [hist; F sum(F,2)];
    dF = hist(end,7)-hist(end-1,7);
    if (dispopt>=3)&cFbool
      for hh = 1:size(hFpart);
	set(hFpart(hh),'XData',1:(size(hist,1)-1),'YData',hist(2:end,hh));
      end
      set(hFdiff,'XData',1:(size(hist,1)-1),'YData',log(diff(hist(1:end,end))));
      drawnow
    end
    if dispopt>=2
      fprintf('\nit:%4i F: B=%.3f A|B=%.3f rho=%.3f D|rho=%.3f C|rho,D=%.3f Y=%.3f F=%.3f dF=%.3f ',it,hist(end,:),dF);
      if dF < 0 % it *should* monotonically increase
	fprintf('\n\n!!!!!\n\n F violation \n\n!!!!!\n');
	alpha, beta, gamma, delta, A,B,C,D
	fprintf('BREAKING (paused)'); pause
	break
      end
    end
  end
  
  % Terminate VBEM iteration loop if number of iterations has been reached, or if the lower bound plateaus.
  if ~(  (it<its) & (dF>1e-14)  )
    break
  end
  
  % Display hint of progress
  if (round(it/50)==it/50) & (dispopt>=1) & ~cFbool, fprintf('.'); end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HYPERPARAMETER UPDATES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if it>=dynhypit
    % ARD for the dynamics
    alpha = k./diag( k*SigA + muA'*muA + 2*SigA*(SA + diag(alpha)*muA' + GA*B')*muA ...
                        + SigA*( (SA+diag(alpha)*muA')*(SA+diag(alpha)*muA')' ...
                                  - 2*GA*B'*( SA+diag(alpha)*muA' )' + GA*( k*SigB+B'*B )*GA' )*SigA ); 
       
    beta  = k./diag( k*SigB + B'*B - 2*B'*muB + muB'*muB); 
  end
  if it>=outhypit
    % ARD for the output process
  
    gamma = p./diag( p*SigC + 2*SigC*MD*D'*rho*muC + muC'*rho*muC...
                      + SigC*( (SC+diag(gamma)*muC')*rho*(SC+diag(gamma)*muC')' ...
                                     - 2*(SC+diag(gamma)*muC')*rho*(D*MD'+muC) + p*MD*SigD*MD' + MD*D'*rho*D*MD')*SigC );
    
    delta = p./diag( p*SigD + D'*rho*D -2*D'*rho*muD + muD'*rho*muD );
  end

  if it>=etchypit
    % Hyperparameters for output noise precision (please see report)
    [pa,pb] = fixedpointsolver(pa,pb,1/p*sum(lnrho,1),1/p*sum(diag(rho),1));
    % Hyperparameters for prior mean and covariance of auxiliary state x_0m (please see report)
    X0m_p = 1/n*sum(cat(2,X0mn{:}),2);
    tmp = zeros(k,k);
    for i=1:n
      tmp = tmp+inv(X0cin{i}) + (X0mn{i}-X0m_p)*(X0mn{i}-X0m_p)';
    end
    X0ci_p = diag(n./diag(tmp));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE HISTORIES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  alphahist = [alphahist alpha];
  betahist = [betahist beta];
  gammahist = [gammahist gamma];
  deltahist = [deltahist delta];

  pahist = [pahist pa];
  pbhist = [pbhist pb];
  X0m_phist = [X0m_phist X0m_p];
  X0ci_phist = [X0ci_phist X0ci_p];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DISPLAY HYP-PROGRESS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if dispopt>=3
    for kk = 1:k
      set(halpha(kk),'XData',1:size(alphahist,2),'YData',-log(alphahist(kk,:)));
      set(hgamma(kk),'XData',1:size(gammahist,2),'YData',-log(gammahist(kk,:)));
    end
    for pinpp = 1:pinp
      set(hbeta(pinpp),'XData',1:size(betahist,2),'YData',-log(betahist(pinpp,:)));
      set(hdelta(pinpp),'XData',1:size(deltahist,2),'YData',-log(deltahist(pinpp,:)));
    end
      set(hrhomean,'XData',1:size(pahist,2),'YData',pahist./pbhist);
      set(hrhovar,'XData',1:size(pahist,2),'YData',log(pahist./(pbhist.^2)));
  end

  
end
if dispopt>=1; fprintf('\n--finished--\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN LEARNT NETWORK FOR ANALYSIS & PERUSAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net = struct('type','Variational State-Space Model (vSSM) with inputs');
net.hyp.pa	= pa;
net.hyp.pb	= pb;
net.hyp.alpha	= alpha;
net.hyp.beta	= beta;
net.hyp.gamma	= gamma;
net.hyp.delta	= delta;
net.hyp.X0m_p	= X0m_p;
net.hyp.X0ci_p	= X0ci_p;
net.param.SigA	= SigA;
net.param.SigB	= SigB;
net.param.SigC	= SigC;
net.param.SigD	= SigD;
net.param.Bb	= Bb;
net.param.Db	= Db;
net.param.SA	= SA;
net.param.SC	= SC;
net.param.GA	= GA;
net.param.MB	= MB;
net.param.MD	= MD;
net.param.G	= G;
net.param.qrhoa	= pa+sum(Tn,2)/2;
net.param.qrhob	= pb+G/2;
net.exp.A	= A;
net.exp.B	= B;
net.exp.AA	= AA;
net.exp.AB	= AB;
net.exp.BB	= BB;
net.exp.rho	= rho;
net.exp.lnrho	= lnrho;
net.exp.C	= C;
net.exp.D	= D;
net.exp.CrhoC	= CrhoC;
net.exp.rhoC	= rhoC;
net.exp.CrhoD	= CrhoD;
net.exp.rhoD	= rhoD;
net.exp.DrhoD	= DrhoD;
net.hidden.X0mn	= X0mn;
net.hidden.Xmn	= Xmn;
net.hidden.X0cin= X0cin;
net.hidden.Xcin	= Xcin;
net.hidden.Ups0n= Ups0n;
net.hidden.Upsn	= Upsn;
net.hist.F	= hist;
net.hist.pa	= pahist;
net.hist.pb	= pbhist;
net.hist.alpha	= alphahist;
net.hist.beta	= betahist;
net.hist.gamma	= gammahist;
net.hist.delta	= deltahist;
net.hist.X0m_p	= X0m_phist;
net.hist.X0ci_p	= X0ci_phist;
