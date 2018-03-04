function [roc_base,VB_base_inter] = get_base_roc(CBDZ_base,info)

stdrange = info.stdrange;
true_inter = info.true_inter;
kk = info.kk;
samples = info.samples;
reps = info.reps;
seed = info.seed;

roc_base = [];
stdind = 0;
for std = stdrange
    stdind = stdind+1;
    CBDZsig_base = CBDZ_base(:,2:end).*(abs(CBDZ_base(:,2:end))>std);
    CBDZsig_base = CBDZsig_base';
    VB_base_inter = CBDZsig_base;
    [se_base sp_base TP_base TN_base FP_base FN_base] = arcroc(true_inter,VB_base_inter);
    roc_base = [roc_base; samples,reps,kk,seed,std,1-sp_base,se_base];
end