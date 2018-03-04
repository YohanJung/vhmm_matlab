%This calculates the second derivative of the log Gamma function.
%
%A hack written with Carl Rasmussen.
%
%Matthew J. Beal GCNU 06/02/01

function res=trigamma(input);

krange = [1:7]';
coef= [-1/12 1/120 -1/252 1/240 -1/132 691/32760 -1/12];

for item = 1:size(input,2);
x = input(item);

if x<=0
   disp(['Error: non-positive argument ' num2str(x) ' in trigamma']);
else if x>=6
   logx=log(x);
   res(item) = 1/x +.5/x^2 -2/x*coef.*krange'*exp(-2*krange*logx);
 else
   res(item) = trigamma(x+1)+1/x^2;
end, end

end