% Calculates the mean and variance of each element in 
% the block-matrix [ A B ; C D ].
% Also returns the significances of each entry.
%
%  [ABCDmean,ABCDvar,ABCDZ] = ABCDioZ(net);
%
% MJB GNCU 11/08/03

function [ABCDmean,ABCDvar,ABCDZ] = ABCDioZ(net);

% Preliminaries
[p k] = size(net.exp.C);
[pinp] = size(net.exp.D,2);

GA = net.param.GA;
  SigA = net.param.SigA;
  SigB = net.param.SigB;
  SigAh = SigA+SigA*GA*SigB*GA'*SigA;
GC = net.param.MD;
  SigC = net.param.SigC;
  SigD = net.param.SigD;
  SigCh = SigC+SigC*GC*SigD*GC'*SigC;
Erhoinv = net.param.qrhob/(net.param.qrhoa-1);

Amean = net.exp.A;
Bmean = net.exp.B;
Cmean = net.exp.C;
Dmean = net.exp.D;

Avar = repmat(diag(SigA)',[k 1]);
Bvar = repmat(diag(SigB)',[k 1]);
Cvar = Erhoinv*diag(SigCh)';
Dvar = Erhoinv*diag(SigD)';

ABCDmean = [Amean Bmean; Cmean Dmean];
ABCDvar = [Avar Bvar; Cvar Dvar];
ABCDZ = ABCDmean./sqrt(ABCDvar);
