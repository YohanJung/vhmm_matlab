function [roc_advc] = get_advc_roc(CBDZ_advc,ind,mu,info)

stdrange = info.stdrange;
true_inter = info.true_inter;
kk = info.kk;
samples = info.samples;
reps = info.reps;
seed = info.seed;

roc_advc = [];
stdind = 0;
for std = stdrange
    stdind = stdind+1;
    CBDZsig_advc = CBDZ_advc(:,2:end).*(abs(CBDZ_advc(:,2:end))>std);
    VB_advc_inter = CBDZsig_advc';
    [se_advc sp_advc TP_advc TN_advc FP_advc FN_advc] = arcroc(true_inter,VB_advc_inter);
    roc_advc = [roc_advc;samples,reps,kk,seed,ind,mu,std,1-sp_advc,se_advc];
end
