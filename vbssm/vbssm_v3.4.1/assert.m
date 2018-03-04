% from Sam 18/05/00
function [] = assert(condition,message)

if nargin == 1,message = '';end
if(~condition) 
  ddd = dbstack;
  if(length(ddd)>1)
    dname=ddd(2).name; 
    dlinestr=['on line ' num2str(ddd(2).line)];
  else
      dname='command line'; 
      dlinestr='';
  end
  fprintf(1,'!!! assert failure (%s)\n    in function %s %s\n',...
      message,dname,dlinestr); 
end
