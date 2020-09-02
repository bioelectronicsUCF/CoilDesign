function Rprox_tot = Res_prox_only3(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)

  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  
  omega = 2 * pi * f;
  
  for i = 1:n
    length_turn(i) = (4*(Dout-w)-8*(i-1)*(w+s));
  end
  
  k = (1+1i)/delta;
  
  Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );

  Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
  Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
%{  
  for i = 1:n
    Hin_(i)=_Hz_sum(s, n, Dout, Din, i, -w/2,1);
    Hout_(i)=_Hz_sum(s, n, Dout, Din, i, w/2,1);
  end
  H=(Hin_+Hout_)/2;
  delH_=Hin_-Hout_  
%}
  H = Hout - (linspace(0,(n-1)/n,n) + 1/2/n) * (Hout-Hin);
  Hin_ = Hout - linspace(1/n,1,n) * (Hout-Hin);
  Hout_ = Hout - linspace(0,(n-1)/n,n) * (Hout-Hin);
  delH = (Hin-Hout)/n;
  delH_ = linspace(1,1,n) * delH;
   

  [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t);
  
  delH_inc = linspace(0,0,n);
  Hin_inc = linspace(0,0,n);
  Hout_inc = linspace(0,0,n);
  
%{  
  for i = 1:n
    for j = 1:n
      if (j==i)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(t+delta*2)-1/(w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(s+delta/2)-2/(s+w-delta/2)+1/(s+2*w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(w/8+s)-1/(s+w*7/8)-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(s+delta/2)-2/(s+w-delta/2)+1/(s+2*w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(w/8+s)-1/(s+w*7/8)-1/(s+w*9/8)+1/(s+15/8*w));
      endif
    end
  end
  
  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/4/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/4/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      endif
    end
  end
%}

  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j) * (1/(t*2+delta*2)-1/(w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j) * 2/pi *(1/(t+w/2)-1/(w-w/4));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j) * (1/(t*2+delta*2)-1/(w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j) * 2/pi *(1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)*2/pi * (1/(s+delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)*2/pi * (1/(w/8+s));
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)*2/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)*2/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)*2/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)*2/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)*2/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)*2/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      end
    end
  end

  w_ = ((Dout - Din)/2)/(2*n-1);
  s_ = w_;
  if s>s_
  else
%  delH_ = delH_ + delH_inc;
  Hin_ = Hin_+Hin_inc;
  Hout_ = Hout_+Hout_inc;
  end
  delH_ = Hin_ - Hout_;


  HE_ = (1-delH_*t)/2/w;  
%  HE = 1*(1-delH*t)/2/w;
%  HE_ = linspace(1,1,n) * HE;
  
  HE_inc = linspace(0,0,n);

  HE_ = HE_ + HE_inc;
  
  alpha = HE_./delH_;
  beta = Hin./delH_;
  
  Rprox_tot=0;
  Rdc_tot=0;
  for i = 1:n
    cothkw=coth(k*w);
    if(isnan(cothkw))
      if(real(k*w)>0)
        cothkw=1;
      else
        cothkw=-1;
      end
    end
    cothtk=coth(t*k/2);
    if(isnan(cothtk))
      if(real(t*k/2)>0)
        cothtk=1;
      else
        cothtk=-1;
      end
    end
    tanhkw=tanh(k*w/2);
    if(isnan(tanhkw))
      if(real(k*w)>0)
        tanhkw=1;
      else
        tanhkw=-1;
      end
    end
    Rprox_per_dc(i) = abs(real(w*t/(2*alpha(i)*w+t)^2*(t*k*cothkw+2*k*w*alpha(i)^2*cothtk+4*alpha(i)+2*t*k*tanhkw*abs(Hin_(i))*abs(Hout_(i))/abs(delH_(i))^2)));%*(beta(i)+1-i)*(beta(i)-i))));%
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  end