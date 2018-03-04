function auc = aucroc(csp,se);
% Area Under ROC curve.
%
% auc = aucroc(csp,se)
%
% Takes complementary specificity and sensitivities.

if ~exist('csp')
  csp = [0 1]';
end
if ~exist('se')
  se = [0 1]';
end


csp = csp(:); % FALSE NEGATIVE RATE
se  = se(:);  % TRUE POSITIVE RATE

pairs = [csp se];

pairs = sortrows(pairs,[1 2]);

if ~all(pairs(1,:)==[0 0])
  %fprintf('Adding 0,0 entry to ROC data\n');
  pairs = [0 0; pairs];
end
if ~all(pairs(end,:)==[1 1])
  %fprintf('Adding 1,1 entry to ROC data\n');
  pairs = [pairs; 1 1];
end

widths = diff(pairs(:,1),1,1);
heights = diff(pairs(:,2),1,1)/2+pairs(1:(end-1),2);

auc = widths'*heights;