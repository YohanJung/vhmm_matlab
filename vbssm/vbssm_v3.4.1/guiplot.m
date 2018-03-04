% function called by guilaunch and guiclickstep to update the graph.
%
% MJB T.O. 22/08/03

function [] = guiplot(nodepos,accept,donor,sigs,labels)

numlinks = size(accept,1);
uniqq = unique(union(accept,donor));
os = .00; rr = sqrt(sum( (max(nodepos,[],1)-min(nodepos,[],1)).^2 ));

colordef none; cla; set(gcf,'doublebuffer','on'); hold on;
for j = 1:numlinks;
  if donor(j)==accept(j)
    if sigs(j)>0 % green for positive self-influence
      plot(nodepos(donor(j),1),nodepos(donor(j),2),'o','markersize',30,'markeredgecolor',[0 1 0]);
    else % red for negative self-influence
      plot(nodepos(donor(j),1),nodepos(donor(j),2),'o','markersize',30,'markeredgecolor',[1 0 0]);
    end
  else
    xd = nodepos(donor(j),1); yd = nodepos(donor(j),2);
    xa  = nodepos(accept(j),1); ya = nodepos(accept(j),2);
    h1 = line( [xd xa], [yd ya] );
    vec = [xa-xd ya-yd]; vec = vec/2;%vec = vec/sqrt(sum(vec.^2)); vec=vec*.05;
    h2 = line( [xd xd+vec(1)] , [yd yd+vec(2)] );
    if sigs(j)>0 % green for positive influence
      set(h1,'color',[0 1 0],'linestyle','-');
      set(h2,'color',[0 1 0],'linestyle','-','linewidth',3);
    else % red for negative influence
      % if you want dotted lines to represent -ve influence, then
      % uncommented next two lines, and comment out further 2 lines.
      %set(h1,'color',[1 0 0],'linestyle','--');
      %set(h2,'color',[1 0 0],'linestyle','--','linewidth',3);
      set(h1,'color',[1 0 0],'linestyle','-');
      set(h2,'color',[1 0 0],'linestyle','-','linewidth',3);
    end
  end
  
end

%plot lovely white blobs;
plot(nodepos(uniqq,1),nodepos(uniqq,2),'o','markersize',20,'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0]);

os=.004;
for j = uniqq';
  %r = text(nodepos(j,1)-os,nodepos(j,2),num2str(j));
  r = text(nodepos(j,1)-os,nodepos(j,2),labels{j});
  set(r,'fontsize',8,'color',[0 0 1]);
end
set(gca,'xtick',[],'ytick',[],'xcolor',[.6 .6 .6],'ycolor',[.6 .6 .6],'position',[.05 .05 .9 .9],'ylim',[-.1 1.1],'xlim',[-.1 1.1])
box(gca,'on')
hold off

