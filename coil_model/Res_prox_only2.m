
function Rprox_tot = Res_prox_only(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)

  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  
  omega = 2 * pi * f;
  
  for i = 1:n
    length_turn(i) = (4*(Dout-w)-8*(i-1)*(w+s));
  end
  
  k = (1+j)/delta;
  
  Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );
  Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
  Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
  H = Hout - (linspace(0,(n-1)/n,n) - 1/2/n) * (Hout-Hin);
  Hin_ = Hout - linspace(1/n,1,n) * (Hout-Hin);
  Hout_ = Hout - linspace(0,(n-1)/n,n) * (Hout-Hin);
  
  delH = (Hin-Hout)/n;
  delH_ = linspace(1,1,n) * delH;


  [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t);
  delH_inc = linspace(0,0,n);
  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(t+delta*2)-1/(w-delta/2));
        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(t+w/2)-1/(w-w/4));

      elseif (abs(i-j)==1)
        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(s+delta/2)-2/(s+w-delta/2)+1/(s+2*w-delta/2));
        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(w/8+s)-1/(s+w*7/8)-1/(s+w*9/8)+1/(s+15/8*w));
      endif
    end
  end
  
  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        Hin_(i) = Hin_(i) + I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hin_(i) = Hin_(i) + I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
        Hout_(i) = Hout_(i) - I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hout_(i) = Hout_(i) - I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
        Hout_(i) = Hout_(i) + I_eddy_edge(j)/4/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hout_(i) = Hout_(i) + I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hin_(i) = Hin_(i) - I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hin_(i) = Hin_(i) - I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
        Hin_(i) = Hin_(i) + I_eddy_edge(j)/4/pi * (1/(s+delta/2)-2/(s+w-delta/2));
        Hin_(i) = Hin_(i) + I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hout_(i) = Hout_(i) - I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hout_(i) = Hout_(i) - I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      endif
    end
  end

  delH_ = delH_ + delH_inc;

  delH_ = Hin_ - Hout_;
  
  HE = 1*(1-delH*t)/2/w;
  HE_ = linspace(1,1,n) * HE;
  
  HE_inc = linspace(0,0,n);

  HE_ = HE_ + HE_inc;
  
  alpha = HE_./delH_;
  beta = Hin./delH_;
  
  Rprox_tot=0;
  Rdc_tot=0;
  for i = 1:n
    Rprox_per_dc(i) = abs(real(w*t/(2*alpha(i)*w+t)^2*(t*k*coth(k*w)+2*k*w*alpha(i)^2*coth(t*k/2)+4*alpha(i)+2*t*k*tanh(k*w/2)*(beta(i)+1-i)*(beta(i)-i))));%*Hin_(i)*Hout_(i)/delH_(i)^2)));
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  
endfunction