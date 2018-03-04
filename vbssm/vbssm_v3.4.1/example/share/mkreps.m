function logynsd = mkreps(y,t,reps,noise)

% Simple script to make a close approximation to the data set presented in
% Zak et al.'s Genome Res. article

% This script assumes switchQ.m is modified to include the line:
%   ts = [540; 60; 10000; 0; 0; 0; 0; 0; 0; 0];
% in the relevant place.
%
% Input:    
%           reps    --- Number of replicate sequences
%           noise   --- 
% Output:   
%           logynsd --- "loged" results
%
% Juan Li 05/30/06

% At each time point, divide the concentrations of each value by its %
% concentration at time 0 (i.e. divide every row by the first row).  % And
% then take the logarithm of the result.  Ignore divide-by-zeros for some
% dimensions.

logyn = log(y)-log(y(ones(500,1),:));

sd = max(logyn,[],1)*noise; % size of noise is noise*max value of each column

logynsd = zeros(500,55,reps);
for rep = 1:reps
  logynsd(:,:,rep) = logyn + randn(size(logyn))*diag(sd);
end

