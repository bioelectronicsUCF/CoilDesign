
function R = Res_tot6(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)
  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  R_coil = Sheet_R_coil * length / w;
  Sheet_R_con;
  R_connection = Sheet_R_con * l_con / w;
  
  rho = Sheet_R_coil*t;
  t_con = rho / Sheet_R_con;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  omega = 2 * pi * f;
  
  w_ = ((Dout - Din)/2)/(2*n-1);
  s_ = w_;
  
  omega_crit = 0;
  %omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil;
  
  if w > t
    if t > 2*delta
      Rskin_tot = R_coil * w * t / (w*t - (w-2*delta)*(t-2*delta));
    else
      Rskin_tot = R_coil * t/delta /(1-exp(-t/delta));
    end
  else
    if w > 2*delta
      Rskin_tot = R_coil * w * t / (w*t - (w-2*delta)*(t-2*delta));
    else
      Rskin_tot = R_coil * w/delta /(1-exp(-w/delta));
    end
  end

  if w > t_con
    if t_con > 2*delta
      Rskin_tot = Rskin_tot + R_connection * w * t_con / (w*t_con - (w-2*delta)*(t_con-2*delta));
    else
      Rskin_tot = Rskin_tot + R_connection * t_con/delta /(1-exp(-t_con/delta));
    end
  else
    if w > 2*delta
      Rskin_tot = Rskin_tot + R_connection * w * t_con / (w*t_con - (w-2*delta)*(t_con-2*delta));
    else
      Rskin_tot = Rskin_tot + R_connection * w/delta /(1-exp(-w/delta));
    end
  end
  
  
  RDC = R_coil + R_connection;
  Rskin_inc = Rskin_tot - RDC;
  Rprox_tot = Res_prox_only15(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f);
  if (omega > omega_crit)
    R = RDC + Rskin_inc + Rprox_tot;
  else
    R = RDC + Rskin_inc + Rprox_tot ;
  end 
  

end
