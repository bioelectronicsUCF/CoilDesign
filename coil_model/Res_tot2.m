function R = Res_tot2(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)
  Sheet_R_coil;
  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  R_coil = Sheet_R_coil * length / w;
  Sheet_R_con;
  R_connection = 2 * Sheet_R_con * l_con / w + Sheet_R_con*(Dout-Din)/w;
  RDC = R_coil + R_connection ;
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  omega = 2 * pi * f;
  
  w_ = ((Dout - Din)/2)/(2*n-1);
  s_ = w_;
  omega_crit = 15 / mu * (w+s) / w^2 * Sheet_R_coil;
  
  if t > 2*delta
    Rskin_tot = RDC * w * t / (w*t - (w-2*delta)*(t-2*delta));
  else
    Rskin_tot = RDC * t/delta /(1-exp(-t/delta));
  end

  Rskin_inc = Rskin_tot - RDC;

  Rprox_tot = Res_prox_only3(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f);
  if (omega > omega_crit)
    R = RDC + Rskin_inc + Rprox_tot;
  else
    R = RDC + Rskin_inc;
  end 

end