function [area_advc] = get_advc_auc(roc_advc)
csp_advc = [];                   se_advc = [];
csp_advc = roc_advc(:,8);        se_advc = roc_advc(:,9);
area_advc = 0;
area_advc = aucroc(csp_advc,se_advc);