function [] = resetr(seed);
if nargin==0
  seed = 0
end
% resets the state of the random number generator matlab 5 style
randn('state',seed);
rand('state',seed);
