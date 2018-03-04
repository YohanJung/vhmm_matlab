data={
  % abc sequences
  'abcabcabcabcabcabcabcabcabcabcabcabc'
  'bcabcabcabcabcabcabcabcabcabcabcabc'
  'bcabcabcabcabcabcabcabcabcabcabcabc'
  'cabcabcabcabcabcabcabcabcabcabcabc'
  'abcabcabc'
  'bcabcabc'
  'abcabcabcabc'
  'bcabcabcabcabcab' 
  'abcabcabcabcabcabc'
  % acb sequences
  'acbacbacbacbacbacbacbacb'
  'bacbacbacbacbacbac'
  'acbacbacbabcbacbacbacbacbacbacbacbacbac'
  'acbacbacbabcbacbacbacbacbacbacbacbacbac'
  'cbacbabcbacbacbacbacbacbacbacbacbac'
  'acbabcbacbacbacbacbacbacbacbacbac'
  'acbacbacb'
  % random sequences
  'aabbbabbbababbbaabbbbababbbaaa'
  'aabbababbababbbbabbaaababaabaa'
  'abbabbbbbbbaabbabbaaaaaabababa'
  'baabaabbabaaaabbabaaabbaabbbaa'
  'bbaaabbababaababbbbbaaabaaabba'
  };

disp('Showing data...');
data
disp(['The data is made of abc sequences, acb sequences, and random ab' ...
      ' sequences']);
fprintf('\n<paused>\n'); pause;
%%
net = vbhmm(data,'abc',13,100,1e-6);
Fv = vbhmm_cF(data,'abc',net);
%%
fig_size = [100,10,1200,800];
fig = figure()
set(fig,'position',fig_size);

subplot(221);
hinton(net.Wa);
s=title('Dirichlet parameters for Transition Matrix');

subplot(223);
hinton(net.Wb);
s=title('Dirichlet parameters for Emission Matrix');

subplot(222);
plot(net.F(:,end));
grid;
ylabel('F');
title('Learning Curve');

subplot(224);
semilogy(diff(net.F(:,end),1));
ylabel('\Delta F');
grid;
xlabel('Iterations of Variational EM');
