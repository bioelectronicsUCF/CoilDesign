% Given a spacing, coil diameters, integer of complete turns, sheet resistance,
% and interconnection length; return effective resistance
function R = Res(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)
  % Calculate the width from the given (integer) number of turns:
  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n+1)^2*(s+w);    
  R_coil = Sheet_R_coil * length / w;
  R_connection = (Sheet_R_coil+Sheet_R_con) * l_con / w + Sheet_R_con*(Dout-Din)/w;
  RDC = R_coil + R_connection ;
  
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
%{
  omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil;
  n_0 = n/4;
  
  Rprox_tot=0; Rdc_tot=0;
  for i = 1:n
    Rprox_per_dc(i) = 0.5 * (2*pi*f/omega_crit)^2 * (i - n_0)^2 / (n - n_0)^2
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
  
  Rprox_tot
  Rdc_tot
  

  
  
 
  omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil;
  f_crit = omega_crit/2/pi
  k = (1+j)/delta;
  
  Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );
  Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
  Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
  delH = (Hin-Hout)/n;
  HE = (1-delH*t)/2/w;
  
  alpha = HE/delH;
  beta = Hin/delH;
  
  Rprox_tot=0;
  Rdc_tot=0;
  for i = 1:n
    Rprox_per_dc(i) = real(w*t/(2*alpha*w+t)^2*(t*k*coth(k*w)+2*k*w*alpha^2+coth(t*k/2)+4*alpha+2*t*k*tan(k*w/2)*(beta+1-i)*(beta-i)));
    R_dc(i) = Sheet_R_coil / w * (4*(Dout-w)-8*(i-1)*(w+s));
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
%}   
  
  Rskin_tot = RDC * t/(delta*(1-e^(-t/delta)));
  Rskin_inc = Rskin_tot - RDC;
  
  R = RDC + Rskin_inc ;

endfunction