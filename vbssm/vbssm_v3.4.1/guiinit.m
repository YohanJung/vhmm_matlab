% Function to initialise the GUI with 'numsigs' significant interactions.
% If a second output is requested, this will randomise the positions.
%
% So to begin with:
%
%   [thresh,nodepos] = guiinit(Z,numsigs)
%
% And subsequently:
%
%             thresh = guiinit(Z,numsigs)
%
% Follow both these up with the 'guilaunch' script.
%
% MJB 21/08/03 T.O.

function [thresh,nodepos] = guiinit(Z,numsigs)

% use numsigs = Inf to obtain all interactions

nodepos = rand(size(Z,1),2);

numsigs = min(numsigs,prod(size(Z)));

if numsigs == 0
  thresh = Inf;
else
  % arrange for 'numsigs' of maximum significance, whether positive or negative.
  sortedZ = flipud(sort(abs(Z(:)))); top10sigs = sortedZ(1:10)
  thresh = sortedZ(numsigs)-1e-10;
end

numsigs = sum(sum(abs(Z)>thresh))
active_genes = sum( (sum(abs(Z)>thresh,1)>0)' | (sum(abs(Z)>thresh,2)>0) )

  