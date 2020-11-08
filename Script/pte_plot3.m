
function pte_plot3(x,pte,xlabel_,ylabel_,xmin,xmax,x_res,ymin,ymax,y_res)
  grid on;
  box off
  fs=9;lw=0.5;
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  plot(x,pte,'ko-','linewidth',lw,'markersize',3,'markerfacecolor','none')
 
  if(abs(1-exist('xmin'))>0)
    xmin=0;
  end
  if(abs(1-exist('xmax'))>0 || xmax=='a')
    xmax=max(xlim);
    xlim([xmin xmax])
  else
    xlim([xmin xmax])
  end

  if(abs(1-exist('ymin'))>0)
    if min(k) > 0
      ymin=0;
    else
      ymin=min(k);
    end
  elseif(ymin=='a')
    ymin=min(k);
  end
  if(abs(1-exist('ymax'))>0 || ymax=='a')
    ymax=max(ylim);
    ylim([ymin ymax])  
  else
    ylim([ymin ymax])  
  end

  if(abs(1-exist('x_res'))>0) || x_res=='a'
    x_res=5;
  else
    set(gca,'xtick',[xmin:(xmax-xmin)/x_res:xmax]);    
  end 
  if(abs(1-exist('y_res'))>0 || y_res=='a')
    y_res=3;
  else
    set(gca,'ytick',[ymin:(ymax-ymin)/y_res:ymax]);    
  end
  
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  aa=get(gca,'position');
  set(gca,'ticklength',[0.001*1/aa(3) 0]); 
  set(gca,'tickdir','out');
  x = get(gca,'position');
  x(3)=x(3);
  set(gca,'position',x);
  
  set(gca,'xminorgrid','off','yminorgrid','off','minorgridlinestyle','-');
  major_gc_=0;
  minor_gc_=0.75;
  majorgridcolor_ = [major_gc_ major_gc_ major_gc_];
  minorgridcolor_ = [minor_gc_ minor_gc_ minor_gc_];
  set(gca,'gridcolor', majorgridcolor_, 'minorgridcolor',minorgridcolor_);
  xlabel(xlabel_);ylabel(ylabel_);
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');  
  grid on;
end

