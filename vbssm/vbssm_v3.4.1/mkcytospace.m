% Function to make files for cytospace format from 
% adjacency matrix and a threshold;
%
% Writes 3 tab-delimited files called 'filename{1,2,3}.cyto'
%
% If you don't specify a threshold, it takes all possible 
% interactions.  A good threshold is sig=3
%
% MJB T.O. 22/08/03

function [data1,data2,data3] = mkcytospace(Z,filename,thresh)

if nargin<3
  thresh=-1;
end

[i j v] = find(Z.*(abs(Z)>thresh));

% i is row index, i.e. acceptor
% j is column index, i.e. donor
% v is significance, may be -ve or +ve

data1 = [j i];

data2 = [j i v>0];

data3 = [j i abs(v)];

fprintf('Writing 3 files: %s{1,2,3}.cyto ...',filename);

dlmwrite([filename '1.cyto'],data1,'\t');
dlmwrite([filename '2.cyto'],data2,'\t');
dlmwrite([filename '3.cyto'],data3,'\t');

fprintf('Done.\n');