function [success] = main(seed)
%%***************************************************
% Show AUC curves with/without ARD prior incorporated, w.r.t varied k

% Support Both Local and Distributed Computing Approaach

% Output: success (0 or 1), indicating function is finished

%%********************************************************
% HOW TO RUN:
% 1. Local Usage:
% >> main(1)
% 
% 2. Remote Usage:(Cluster)
%  >> matlab    % open matlab program
%     
%  >> pppp = {'/home/csgrad/juanli/work/toolkits/vbssm_v3.4.1',...
%           '/home/csgrad/juanli/work/code'};
%     % pppp is your working path and has to be a cell!
%  
%  >> s = dfeval(@main,num2cell(1),'PathDependencies',pppp)
%     % Evaluate function using cluster, with
%     % seed provided in the cell arrays.
%%*************************************************************

local = 1;
reprange = [1];
arc_ind = 5%  the arc of interest, i.e. a prior will be added on this arc;
kkrange = [1:16];
success = 0;
resetr(seed);
if local % Local setting
    addpath('c:\matlab701\work\vbssm_v3.4.1');
    addpath(genpath('C:\MATLAB701\work\0722\code'));

    load('true_inter.mat');
    load('y.mat');
    load('t.mat');

else % Remote Setting
    %     Load ODE data and true interaction

    addpath('/home/csgrad/juanli/work/toolkits/vbssm_v3.4.1');
    addpath( genpath('/home/csgrad/juanli/work/code') );

    load('/home/csgrad/juanli/work/code/share/true_inter.mat');
    load('/home/csgrad/juanli/work/code/share/y.mat');
    load('/home/csgrad/juanli/work/code/share/t.mat');
end

%%%%%%%%%%%%%% Global Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdrange = [0.5:0.5:10];
mu = 2; samples = 12;reps = 1;
noise = 0.1;  % set noise to 10%
maxreps = 200; % set number of replicates
% ideal time window 600-3000(min) == 10-50(hr)
time_window = [600 3480];
vbits = 300;

info.stdrange =stdrange ;
info.true_inter = true_inter;
info.seed = seed;
info.reps = reps;
info.samples = samples;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a close approximation to the data set presented in
% Zak et al.'s Genome Res. article
logynsd = mkreps(y,t,maxreps,noise);
% data sampling, time_window defines start & end time

all_logynsd =  sampling(logynsd,t,time_window,samples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampled_logynsd = all_logynsd(:,:,1:reps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T mRNA reps] = size(sampled_logynsd);
yn = {}; Tn = {}; inpn ={};
for rep = 1:reps
    yn{rep} = sampled_logynsd(2:T,:,rep)';
    Tn = size(yn{rep},2);
    inpn{rep} = [ones(1,Tn);sampled_logynsd(1:(T-1),:,rep)'];
end

p = size(yn{1},1);		% take first obs seq to set the observ. dim
pinp = size(inpn{1},1);		% take first input sequence to set input dim

[g_rows,g_cols] = find_arcs(true_inter);


% Find the single_arc of interest
col = g_cols(arc_ind);   row = g_rows(arc_ind);

auc_base = [];
auc_advc = [];

kkind = 0;
for kk =  kkrange
    info.kk = kk;
    kkind = kkind +1;
    
    muA= zeros(kk,kk); muB = zeros(kk,pinp); muC= zeros(p,kk); muD = zeros(p,pinp);
    ZeroPrior{1} = muA;   ZeroPrior{2} = muB;  ZeroPrior{3} = muC; ZeroPrior{4} = muD;
    %   Train Basic VBSSM, Without PRIOR Info
    VB_base_net = vssminpn(yn,inpn,kk,vbits,0,0,ZeroPrior,[Inf Inf Inf]);

    %to compute significance of these distributions for structure determination.
    [CBDmean_base,CBDvar_base,CBDZ_base,CBDZsig_base] = CBDioZ(VB_base_net);
    %to compute ROC and AUC measurements of an inferred structure from a known structure.
    [roc_base,VB_base_inter] = get_base_roc(CBDZ_base,info);
    [area_base] = get_base_auc(roc_base);
    
    auc_base = [auc_base,area_base ];
    fprintf('kk=%2i: Auc=%2.5f (w/o ARD prior)\n',kk,area_base);

    muA= zeros(kk,kk); muB = zeros(kk,pinp); muC= zeros(p,kk);muD = zeros(p,pinp);
    PriMu{1} = muA;   PriMu{2} = muB;  PriMu{3} = muC;
    tmp2 = zeros(p,p);
    % Add ARD Prior Info
    tmp2(row,col) = sign(true_inter(row,col))*mu;
    muD(:,2:end) = tmp2';
    PriMu{4} = muD;

    % Train Advanced VBSSM, With Prior Info
    VB_advc_net = vssminpn(yn,inpn,kk,vbits,0,0,PriMu,[Inf Inf Inf]);
    % to compute significance of these distributions for structure determination.
    [CBDmean_advc,CBDvar_advc,CBDZ_advc,CBDZsig_advc] = CBDioZ(VB_advc_net);
    %to compute ROC and AUC measurements of an inferred structure from a known structure.
    [roc_advc] = get_advc_roc(CBDZ_advc,arc_ind,mu,info);
    [area_advc] = get_advc_auc(roc_advc);
    auc_advc = [auc_advc,area_advc ];
    fprintf('       Auc=%2.5f (with ARD prior)\n',area_advc);
    %**********************************************
end

figure(1);
plot(auc_base,'.-b');hold on,
plot(auc_advc,'.-r');hold off;
xlabel('K, dimension of hidden space');
ylabel('auc');
set(gca,'XTick',[1:16]);
set(gca,'XTickLabel',{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16'});
    
legend('Auc, w/o ARD prior added','Auc, with ARD prior added');
success = 1;

disp('End');

