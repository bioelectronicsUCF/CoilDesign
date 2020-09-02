function status=kplot_(x,k,xlabel_,xmin,xmax,x_res,ymin,ymax,y_res)
  hold on;
  grid on;
  box off
  
  if(abs(1-exist('xmin'))>0)
    xmin=min(x);
  end
  if(abs(1-exist('xmax'))>0||xmax==0)
    xmax=max(max(x),max(xlim));
  end
  if(abs(1-exist('x_res'))>0)
    x_res=3;
  end
  if(abs(1-exist('ymin'))>0)
    ymin=0;
  end
  if(abs(1-exist('ymax'))>0)
    ymax=max(k);
  end
  if(abs(1-exist('y_res'))>0)
    y_res=3;
  end
  
  fs=9;lw=0.25;
  plot(x,k,'ro-','linewidth',lw,'markersize',3,'markerfacecolor','none')

  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  aa=get(get(gca,'parent'),'position');
  aaa=get(gca,'position')';
  set(gca,'ticklength',[0.005*1/aaa(3) 0]); 
  set(gca,'tickdir','out');
%  set(gca,'xticklabel',{''});
  set(gca,'xtick',[xmin:(xmax-xmin)/x_res:xmax]);
%  set(gca,'yticklabel',{''});
  set(gca,'ytick',[ymin:(ymax-ymin)/y_res:ymax]);
  xlim([xmin xmax]);
  ylim([ymin ymax]);

  

  x = get(gca,'position');
  x(3)=x(3);
  set(gca,'position',x)
  
  set(gca,'xminorgrid','off','yminorgrid','off','minorgridlinestyle','-');
  major_gc_=0;
  minor_gc_=0.75;
  majorgridcolor_ = [major_gc_ major_gc_ major_gc_];
  minorgridcolor_ = [minor_gc_ minor_gc_ minor_gc_];
  set(gca,'gridcolor', majorgridcolor_, 'minorgridcolor',minorgridcolor_);
  xlabel(xlabel_);ylabel("PTE (%)");
%  legend('Cal','Sim');
%{  
  for i = ymin:(ymax-ymin)/y_res:ymax
    line([xmin -xmax/30*x(4)/x(3)*(ymax-ymin)/(xmax-xmin)],[i i],'linewidth',2,'linestyle','-','color','k');
    text(xmin,i,[num2str(i) " "],'horizontalalignment','right')
  end
  line([xmin xmin],[ymin ymax], 'linewidth',2,'linestyle','-','color','k'); 

  text(xmin,-0.1,"0",'horizontalalignment','left','verticalalignment','top','fontsize',fs)
  text(xmax,-0.1,num2str(xmax),'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  
  for i = 0:xmax/3:xmax
    line([i i],[0 max(ylim)/15],'linewidth',2,'linestyle','-','color','k');
  end
  line([0 max(xlim)],[0 0], 'linewidth',2,'linestyle','-','color','k'); 
%}
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  
endfunction
