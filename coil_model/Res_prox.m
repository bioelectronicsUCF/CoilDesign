% Given a spacing, coil diameters, integer of complete turns, sheet resistance,
% and interconnection length; return effective resistance
function R = Res_prox(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)
  % Calculate the width from the given (integer) number of turns:
  Sheet_R_coil;
  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  R_coil = Sheet_R_coil * length / w;
  Sheet_R_con;
  R_connection = 2 * Sheet_R_con * l_con / w + Sheet_R_con*(Dout-Din)/w;
  RDC = R_coil + R_connection ;
  
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5
  
  omega = 2 * pi * f
  omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil
  if omega>omega_crit*5
%    omega = 5 * omega_crit;
  end
  
  for i = 1:n
    length_turn(i) = (4*(Dout-w)-8*(i-1)*(w+s));
  end
  
%  n_0 = n/4;
%{
  Rprox_tot=0; Rdc_tot=0;
  for i = 1:n
    Rprox_per_dc(i) = 0.5 * (omega/omega_crit)^2 * (i - n_0)^2 / (n - n_0)^2;
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  
  Rprox_tot
  Rdc_tot
%}  

%{
  Rprox_tot=0; Rdc_tot=0;
  for i = 1:n
    Rprox_per_dc(i) = 0.14 * w^4 * t * (1-e^(-t/delta))  / delta^3 / (s+w)^2 * (i - n_0)^2 / (n - n_0)^2;
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  
  Rprox_tot
  Rdc_tot
%}

%
  omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil;
  f_crit = omega_crit/2/pi;
  k = (1+j)/delta;
  
%  Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );
  Hc=_Hz_sum(s, n, Dout, Din, n, -Din/2-w/2,1);

  Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
  Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
  H = Hout - (linspace(0,1,n) - 1/2/n) * (Hout-Hin)
  
  delH = (Hin-Hout)/n;
  delH_ = linspace(1,1,n) * delH;

  w,delta,H,Sheet_R_coil,omega,n,t
  [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t)
  

  
  delH_inc = linspace(0,0,n);
  for i = 1:n
    for j = 1:n
      if (abs(i - j)<1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(t+delta*2)-1/(w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(t+w/2)-1/(w-w/4));
        
      elseif (abs(i-j)==1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(abs(i-j-1)*(w+s)+s+delta/2)-2/(abs(i-j-1)*(w+s)+s+w-delta/2)+1/(abs(i-j-1)*(w+s)+s+2*w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(abs(i-j-1)*(w+s)+w/8+s)-1/(abs(i-j-1)*(w+s)+s+w*7/8)-1/(abs(i-j-1)*(w+s)+s+w*9/8)+1/(abs(i-j-1)*(w+s)+s+15/8*w));
         delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * ((s+delta/2)/((s+delta/2)^2+(t/2)^2)-(s+w-delta/2)/((s+w-delta/2)^2+(t/2)^2));
         delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * ((s+delta/2)/((s+delta/2)^2+(t/2)^2)-(s-delta/2+w*7/8)/((s-delta/2+w*7/8)^2+(t/2)^2))

      endif
    end
  end
  delH_ = delH_ + delH_inc;
  
  HE = (1-delH*t)/2/w;
  HE_ = linspace(1,1,n) * HE;
  
  HE_inc = linspace(0,0,n);
%{
  for i = 1:n
    for j = 1:n
      if (abs(i-j)<1)
        HE_inc(i) = HE_inc(i) + I_eddy_edge(j) /pi/w * atan(w/t) + I_eddy_center (j) / (w/4) /2
      elseif (abs(i-j)==1)
        HE_inc(i) = HE_inc(i) - I_eddy_edge(j) /pi/w/2 * (atan((w+s)/(t/2))-atan(s/(t/2)))
      endif
    end
  end
%}
  HE_ = HE_ + HE_inc;
  
  alpha = HE_./delH_
  beta = Hin./delH_
  
  Rprox_tot=0;
  Rdc_tot=0;
  for i = 1:n
    coth_(k*w),coth_(t*k/2),tan(k*w/2), k*w/2;
    Rprox_per_dc(i) = abs(real(w*t/(2*alpha(i)*w+t)^2*(t*k*coth_(k*w)+2*k*w*alpha(i)^2+coth_(t*k/2)+4*alpha(i)+2*t*k*tan_(k*w/2)*(beta(i)-n+i)*(beta(i)-n-1+i))))
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  RDC=Rdc_tot;
  
  Rprox_tot;
  Rprox_tot;
    if w > 2*delta
    if t > 2*delta
      Rskin_tot = RDC * w * t / (w*t - (w-2*delta)*(t-2*delta));
    else
      Rskin_tot = RDC * t/(delta*(1-e^(-t/delta)));
    end
  else
    if t > 2*delta
      Rskin_tot = RDC * w/(delta*(1-e^(-w/delta)));
    else
      Rskin_tot = RDC;
    end
  end  
  Rskin_inc = Rskin_tot - RDC;

  cut_=[1.5, 7, 50];
  co_=[0.5,0.3];
  if s > cut_(3)*w
    R = Rdc_tot + Rskin_inc + 0* Rprox_tot
  elseif s > cut_(2)*w
    R = Rdc_tot + Rskin_inc + co_(2)*(s-cut_(3)*w)/(cut_(2)*w-cut_(3)*w)*Rprox_tot
  elseif s > cut_(1)*w
    R = Rdc_tot + Rskin_inc + (co_(1)*(cut_(2)*w-s)+co_(2)*(s-cut_(1)*w))/(cut_(2)*w-cut_(1)*w)* Rprox_tot
  else
    R = Rdc_tot +Rskin_inc + Rprox_tot
  end
  
  function y=coth_(x)
    if isnan(coth(x))
      y=1;
    else
      y=coth(x);
    end
  endfunction
  
  function y=tan_(x)
    if isnan(tan(x))
      if imag(x)>5
        y=1i;
      elseif imag(x)<5
        y=-1i;
      else
        y=0
      end
    else
      y=tan(x);
    end
  endfunction
    
  
endfunction