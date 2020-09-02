function R_turns = Res_tot_turns(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)

  [Rdc_turns, Rprox_turns] = Res_prox_only2_turns(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f);
  
  w = ((Dout - Din)/2 - s*(n-1))/n;
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  omega = 2 * pi * f;
  omega_crit = 3.1 / mu * (w+s) / w^2 * Sheet_R_coil;
  
  if t > 2*delta
    Rskin_turns = Rdc_turns * w * t / (w*t - (w-2*delta)*(t-2*delta));
  else
    Rskin_turns = Rdc_turns * t/delta /(1-e^(-t/delta));
  end

  if (omega > omega_crit)
    R_turns = Rskin_turns + Rprox_turns;
  else
    R_turns = Rskin_turns;
  end 

endfunction