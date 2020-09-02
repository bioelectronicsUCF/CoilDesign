function Cp_turns = Capp_turns(s, n, od, id, epsilon_r, t)

  w = ((od - id)/2 - s*(n-1))/n;
  
  gamma = 0.577216;
  epsilon_0 = 8.8541878e-18;  % permittivity of vacuum, F/um
  epsilon = epsilon_0 * epsilon_r;  % permittivity of insulator between metals, F/um
  
  for i = 1 : n-1
    len_in = 4*od-8*w-8*(i-1)*(w+s);
    rr = (s+w*2)/s;
    C_parallel_plate = epsilon / s * t * len_in;
    C_fringe = pi * epsilon * len_in /rr * ellipke((rr^2-1)/rr^2) / ( 1/(rr+1)*(2*(gamma+log(4))*(ellipke((rr-1)^2/(rr+1)^2))-pi*(p_deriv(3,[0.5,0.5,1,(rr-1)^2/(rr+1)^2])+p_deriv(2,[0.5,0.5,1,(rr-1)^2/(rr+1)^2]))) + log((rr+1)/(rr-1))/rr*ellipke((rr^2-1)/rr^2) );
    Cp_turns(i) = C_parallel_plate + C_fringe;
  end
  
endfunction