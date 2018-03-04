function  [area_base] = get_base_auc(roc_base)
%************* Compute area Under Basic_ROC curve *************
csp_base = [];                   se_base = [];
csp_base = roc_base(:,6);        se_base = roc_base(:,7);
area_base = 0;
area_base = aucroc(csp_base,se_base);