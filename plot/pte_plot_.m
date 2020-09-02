function state = pte_plot_(x,y,xlabel_,ylabel_,title_,path_);
  hold on;
  grid on;
  
  fs=20;lw=3;
  set(gca,'fontsize',fs,'linewidth',lw);

  plot(x,y,'ko-','linewidth',3,'markersize',15)

  xmax=max(xlim)*1.15;
  xmin=min(xlim);
  xlim([xmin xmax]);
%{
  set(gca,'xticklabel',{''});
  text(xmin,-0.02*max(ylim),"0",'horizontalalignment','left','verticalalignment','top','fontsize',fs)
  text(max(xlim),-0.02*max(ylim),num2str(max(xlim)),'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  set(gca,'xtick',[0:xmax/4:xmax]);
  set(gca,'yticklabel',{''});
  ymax=ceil(max(ylim)*10)/10;
  set(gca,'ytick',[0:ymax/3:ymax]);
  text(xmin,0,[num2str(xmin) " "],'horizontalalignment','right','verticalalignment','bottom','fontsize',fs)
  text(xmin,ymax,[num2str(ymax) " "],'horizontalalignment','right','verticalalignment','top','fontsize',fs)
  ylim([0 ymax])

  x = get(gca,'position');
  set(gca,'xminorgrid','off','yminorgrid','off','minorgridlinestyle','-');
  major_gc_=0;
  minor_gc_=0.7;
  grid on;
  majorgridcolor_ = [major_gc_ major_gc_ major_gc_];
  minorgridcolor_ = [minor_gc_ minor_gc_ minor_gc_];
  set(gca,'gridcolor', majorgridcolor_, 'minorgridcolor',minorgridcolor_);
  set(gca,'xminortick','off','yminortick','off');
  set(gca,'ticklength',[0.05 0.1]); 
  for i = 0:ymax/15:ymax
    line([xmin max(xlim)/30*x(4)/x(3)],[i i],'linewidth',2,'linestyle','-','color','k');
  end   
  for i = 0:ymax/3:ymax
    line([xmin max(xlim)/15*x(4)/x(3)],[i i],'linewidth',2,'linestyle','-','color','k');
  end
  line([xmin xmin],[0 max(ylim)], 'linewidth',2,'linestyle','-','color','k'); 
  
  line([xmin max(xlim)],[0 0], 'linewidth',2,'linestyle','-','color','k'); 
%}
%  legend('Calculation', 'Simulation')
%  legend('location','southeast')
  
  xlabel(xlabel_)
  ylabel(ylabel_)
  ylim([0 max(ylim)]);
  saveas(gca,[path_ title_ ".png"])

endfunction
