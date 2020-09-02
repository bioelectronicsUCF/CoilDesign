function status=rodplot_xy_change(sorted_,fs,xmin,xmax,zmin,zmax,xlabel_,k_,x_,z_,sim_n,zres,xres)
  hold on;grid on;
%{  
  plot(sorted_(1,1:end),sorted_(2,1:end),'ko','markersize',3);
  plot(sorted_(1,1:end),sorted_(3,1:end),'ko','markersize',3);
  plot(sorted_(1,1:end),sorted_(4,1:end),'ko','markersize',3);

  semilogy(sorted_(2,1:end),1./sorted_(1,1:end),'k');
  semilogy(sorted_(3,1:end),1./sorted_(1,1:end),'k:');
  semilogy(sorted_(4,1:end),1./sorted_(1,1:end),'k:');
%}
  zf =sorted_(1,1):(sorted_(1,end)-sorted_(1,1))/5:sorted_(1,end);
  xf2  =interp1 (sorted_(1,1:end),sorted_(2,1:end),zf,'pchip');
  xf3  =interp1 (sorted_(1,1:end),sorted_(3,1:end),zf,'pchip');
  xf4  =interp1 (sorted_(1,1:end),sorted_(4,1:end),zf,'pchip');
  
  zf_=sorted_(1,1):(sorted_(1,end)-sorted_(1,1))/99:sorted_(1,end);
  xf2_  =interp1 (zf,xf2,zf_,'pchip');
  xf3_  =interp1 (zf,xf3,zf_,'pchip');
  xf4_  =interp1 (zf,xf4,zf_,'pchip');
  
  zf__=sorted_(1,1):(sorted_(1,end)-sorted_(1,1))/99:sorted_(1,end);
  xf2__  =interp1 (zf_,xf2_,zf__,'pchip');
  xf3__  =interp1 (zf_,xf3_,zf__,'pchip');
  xf4__  =interp1 (zf_,xf4_,zf__,'pchip');

  zf_=sorted_(1,1):(sorted_(1,end)-sorted_(1,1))/99:sorted_(1,end);
  xf2_  =interp1 (zf__,xf2__,zf_,'pchip');
  xf3_  =interp1 (zf__,xf3__,zf_,'pchip');
  xf4_  =interp1 (zf__,xf4__,zf_,'pchip');
  
  map_s = [183 36 103;
       246 139 31;
       102 45 145;
       82 79 161;
       0 166 81]/255;
       
  for i = 1:32
    st = floor(i/8.01);
    en = ceil(i/8.01);
    map_n(i,:) = map_s(st+1,:) * (en-i/8) + map_s(en+1,:) * (i/8-st);
  end

  colormap(flipud(map_n))  
  caxis([0 1])
  c_ = ones(1,100);
  patch([zf_,fliplr(zf_),zf_,fliplr(zf_)], [xf3_,fliplr(xf2_),xf2_,fliplr(xf4_)], [c_*0.95,c_*1,c_*1,c_*0.95],'linestyle','none');
  plot(zf_,xf2_,'b','linewidth',2);
  plot(zf_,xf3_,'k--');
  plot(zf_,xf4_,'k--');

shading faceted
 
  if(abs(exist("xmin")-1)>0)
    xmin = min(xlim);
  end
  if(abs(exist("zres")-1)>0)
    zres = 5;
  end
  if(abs(exist("xres")-1)>0)
    xres = 5;
  end
  if(abs(exist("xmax")-1)>0)
    xmax = max(xlim);
  end
  if(abs(exist("zmin")-1)>0)
    zmin = min(zlim);
  end
  if(abs(exist("zmax")-1)>0)
    zmax = max(zlim);
  end
  if(abs(exist("xlabel_")-1)>0)
    xlabel_ = '';
  end
  if(abs(exist("fs")-1)>0)
    fs = 9;
  end  
  lw=0.5;


  set(gca,'fontsize',fs,'linewidth',lw,'fontname','Arial');
  aa=get(get(gca,'parent'),'position');
  aaa=get(gca,'position')'
  set(gca,'ticklength',[0.0025*1/aaa(3) 0])
  set(gca,'tickdir','out');

  set(gca,'ytick',[xmin:(xmax-xmin)/xres:xmax]);
  set(gca,'xtick',[zmin:(zmax-zmin)/zres:zmax]);
  
  ylim([xmin xmax])
  xlim([zmin zmax])
  %ylabel(["Rod = " _format_(z) "z\n"])

%  xticks = get(gca,"ytick");
%  ylabels = arrayfun (@(x) format_ (x), yticks, "uniformoutput",false);
%  set(gca,"xticklabel",ylabels);
  
  ylabel(xlabel_);
  xlabel("d / OD{RX}");
    
  function a = format_(x)
    if x>=10
      a = sprintf("%.0f",x);
    elseif x>=1
      a = sprintf("%.0f",x);
    else
      a = sprintf("%.0f",x);
    endif
  endfunction
  function b = _format_(x)
    if x>1
      b = sprintf("%.0f",1/x);
    elseif x>0
      b = sprintf("%.0f",1/x);
    else
      b = sprintf("%.0f",0);
    endif
  endfunction
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','Arial');

endfunction
