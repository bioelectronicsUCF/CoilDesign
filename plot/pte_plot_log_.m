function state = pte_plot_log(x,y,xlabel_,ylabel_,title_,path_,x_div);

  semilogx(x,y,'ko','linewidth',3,'markersize',15)

  fs=40;lw=3;
  set(gca,'fontsize',fs,'linewidth',lw);
  xmax=max(xlim);
  xlim([min(xlim) xmax]);
  set(gca,'xticklabel',{''});
  text(min(xlim),-0.02*max(ylim),num2str(min(xlim)),'horizontalalignment','left','verticalalignment','top','fontsize',fs)
  text(max(xlim),-0.02*max(ylim),num2str(max(xlim)),'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  set(gca,'xtick',[0:xmax/x_div:xmax]);
  set(gca,'yticklabel',{''});
  ymax=ceil(max(ylim)*10)/10;
  set(gca,'ytick',[0:ymax/3:ymax]);
  text(min(xlim),0,"0 ",'horizontalalignment','right','verticalalignment','bottom','fontsize',fs)
  text(min(xlim),ymax,[num2str(ymax) " "],'horizontalalignment','right','verticalalignment','top','fontsize',fs)
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

%{  
  for i = 0:xmax/15:xmax
    line([i i],[0 max(ylim)/30],'linewidth',2,'linestyle','-','color','k');
  end   
  for i = 0:xmax/3:xmax
    line([i i],[0 max(ylim)/15],'linewidth',2,'linestyle','-','color','k');
  end
%}
  line([0 max(xlim)],[0 0], 'linewidth',2,'linestyle','-','color','k'); 
  set(gca,'xtickmode','auto')
  set(gca,'xticklabel','') 
  
  xlabel(xlabel_)
  ylabel(ylabel_)
  ylim([0 max(ylim)]);
  saveas(gca,[path_ title_ ".png"])

endfunction

