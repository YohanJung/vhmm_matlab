%result = lndet(Sigma);
%
%Returns the natural logarithm of the determinant of a
%positive definite matrix, more accurately than the naive way.
%
%Matthew J. Beal 01/01/01

function [result] = lndet(Sigma);

result = 2*sum(log(diag(chol(Sigma))),1);