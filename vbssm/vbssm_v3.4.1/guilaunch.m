% this returns the row/col/value of all those who absolute > 3;

if ~exist('thresh')
  thresh=3;
end

[accept donor sigs] = find(Z.*(abs(Z)>thresh));
[y inds] = sort(abs(sigs)); inds = fliplr(inds);
sigs = sigs(inds);
accept = accept(inds); 
donor = donor(inds);

numint = size(sigs,1);

if ~exist('nodepos')
  nodepos = rand(size(Z,1),2);
end

%if ( size(nodepos,1) ~= size(Z,1) );
%  disp('mismatch in loaded positions and your Z matrix --- randomising positions');
%  nodepos = rand(size(Z,1),2);
%end


if exist('labelint') & ~exist('labels');
  for nn = 1:size(labelint);
    labels{nn} = labelint(nn);
  end
end

% Only overwrite if you don't have labels in memory.
if ~exist('labels');
  % If you specify a vector of numbers to label the genes with, we'll use them.
  if exist('labelint') & ~exist('labels');
    for nn = 1:size(labelint);
      labels{nn} = labelint(nn);
    end
  else % otherwise we'll default to a new set, just boring numbers
    for nn = 1:size(nodepos,1)
      labels{nn} = num2str(nn);
    end
  end
end

% e.g. you could have labels{1} = 'A', such that the first gene is labelled 'A'.  Useful for hidden state.

guiplot(nodepos,accept,donor,sigs,labels) ;

set(gca,'color','none');
set(gca,'buttondownfcn','[nodepos,closest] = feval(@guiclickstep,nodepos,accept,donor,sigs,labels,0,[]);')
