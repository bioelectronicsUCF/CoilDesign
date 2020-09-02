function state = pte_plot(x,y,xlabel_,ylabel_,title_,path_);

  plot(x,y,'linewidth',5)

  fs=40;lw=3;
  set(gca,'fontsize',fs,'linewidth',lw);
  xmax=max(xlim);
  xlim([0 xmax]);
  set(gca,'xticklabel',{''});
  text(0,-0.02*max(ylim),"0",'horizontalalignment','left','verticalalignment','top','fontsize',fs)
  text(max(xlim),-0.02*max(ylim),num2str(max(xlim)),'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  set(gca,'xtick',[0:xmax/4:xmax]);
  set(gca,'yticklabel',{''});
  ymax=ceil(max(ylim)*10)/10;
  set(gca,'ytick',[0:ymax/3:ymax]);
  text(-0,0,"0 ",'horizontalalignment','right','verticalalignment','bottom','fontsize',fs)
  text(-0,ymax,[num2str(ymax) " "],'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  ylim([0 ymax])

  x = get(gca,'position');
  set(gca,'xminorgrid','on','yminorgrid','on','minorgridlinestyle','-');
  major_gc_=0;
  minor_gc_=0.7;
  grid on;
  majorgridcolor_ = [major_gc_ major_gc_ major_gc_];
  minorgridcolor_ = [minor_gc_ minor_gc_ minor_gc_];
  set(gca,'gridcolor', majorgridcolor_, 'minorgridcolor',minorgridcolor_);
  set(gca,'xminortick','on','yminortick','on');
  set(gca,'ticklength',[0.05 0.1]); 
  for i = 0:ymax/15:ymax
    line([0 max(xlim)/30*x(4)/x(3)],[i i],'linewidth',2,'linestyle','-','color','k');
  end   
  for i = 0:ymax/3:ymax
    line([0 max(xlim)/15*x(4)/x(3)],[i i],'linewidth',2,'linestyle','-','color','k');
  end
  line([0 0],[0 max(ylim)], 'linewidth',2,'linestyle','-','color','k'); 
  
  for i = 0:xmax/20:xmax
    line([i i],[0 max(ylim)/30],'linewidth',2,'linestyle','-','color','k');
  end   
  for i = 0:xmax/4:xmax
    line([i i],[0 max(ylim)/15],'linewidth',2,'linestyle','-','color','k');
  end
  line([0 max(xlim)],[0 0], 'linewidth',2,'linestyle','-','color','k'); 

  
  xlabel(xlabel_)
  ylabel(ylabel_)
  ylim([0 max(ylim)]);
  saveas(gca,[path_ title_ ".png"])

endfunction
