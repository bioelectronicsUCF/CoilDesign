% Given a spacing, coil diameters, turn count, insulator relative permittivity,
% insulator height, metal thickness; return effective capacitance
function [C, C_substrate, C_mm_h, C_mm_v] = distributedCap2(s, n, od, id, epsilon_r_sub, epsilon_r_mm, epsilon_r_mmf, d, h, t, c_sub_per_area, g_sub_per_area, f)


  % Calculate the width from the given (integer) number of turns:
  w = ((od - id)/2 - s*(n-1))/n;
  
  gamma = 0.577216;
  epsilon_0 = 8.8541878e-18;  % permittivity of vacuum, F/um

  epsilon_sub = epsilon_0 * epsilon_r_sub;  % permittivity of insulator between metal and substrate, F/um
  epsilon_mm = epsilon_0 * epsilon_r_mm;  % permittivity of insulator between metals, F/um
  epsilon_mmf = epsilon_0 * epsilon_r_mmf;
  
  % First, calculate the capacitance between metal and substrate
  % Need to know the area of each turn
  Area(1) = od^2 - (od-2*w)^2 - s*w;  %area of outermost turn
  for i = 2:n
    Area(i) = (od-2*(i-1)*(w+s))^2 - (od-2*(i-1)*(w+s)-2*w)^2 - s*w;  %metal area of each turn
  end

  % Next is the length of each turns along the center of the metal trace
  Length_c(1) = (od-w/2) + (od-w) + (od-w) + (od-1.5*w-s);
  for i = 2:n
    Length_c(i) = Length_c(1) - 8*(i-1)*(w+s) + s;
  end
  Length_c_tot = 4*n*od - 4*n*w - (2*n-1)^2*(w+s) + w;

  % From the length, calculate the variables for calculating coefficient of capacitance for each turns
  r = Length_c / Length_c_tot;  %ratio of the each length of the turn to total length
  % Do summation of h from 1 turn to i turn and name as s(i)
  for i = 1:n
    sum_r(i) = 0;
    for j = 1:i
      sum_r(i) = sum_r(i) + r(j);
    end
  end

  % not using fringe term which is C_f = epsilon * ( (2*pi) / log(1+2*h/t+sqrt(2*h/t*(2*h/t+2))) - t/2/h );  

  % Now using areas and ds for each turn, we can calculate the capacitance between metal and substrate by adding portions of each turn
  
  b(1) = (1-r(1)/2)^2;
  for i = 2 : n
    b(i) = (1-sum_r(i-1)-r(i)/2)^2;
  end
  C_substrate = 0;
  if g_sub_per_area > 0
    for i = 1 : n
      C_turns_sub(i) = epsilon_sub / h * Area(i);
      C_sub(i) = Area(i) * c_sub_per_area;
      R_sub(i) = 1 / Area(i) / g_sub_per_area;
      C_substrate = C_substrate + b(i) * 1/( 1/C_turns_sub(i)+(2i*pi*f*R_sub(1))/(2i*pi*f*R_sub(1)*C_sub(i)+1) );
    end
  end
    

  C_mm_v = epsilon_mm / d * w^2 * b(1);
  for i = 2 : n
    C_mm_v = C_mm_v + epsilon_mm / d * w^2 * b(i);
  end
    
  % Second, calculate the capacitance between the metals of adjacent turns
  % Need to know the length of each turns along the innermost boundaries of each turns
  Length_c;
  Length_in = Length_c - 4*w;
  Length_out = Length_c + 4*w;

  % Now using lengths and ds for each turn, we can calculate the capacitance between metals of adjacent turns by adding portions of each turn
  C_mm_h = 0;
  for i = 1 : n-1
    for j = i+1  : n
     rr = (s*abs(j-i)+w*2)/(s*abs(j-i));
%      C_mm_h_DC = pi * epsilon_mmf * Length_out(j) /rr * ellipke((rr^2-1)/rr^2) / ( 1/(rr+1)*(2*(gamma+log(4))*(ellipke((rr-1)^2/(rr+1)^2))-pi*(p_deriv(3,[0.5,0.5,1,(rr-1)^2/(rr+1)^2])+p_deriv(2,[0.5,0.5,1,(rr-1)^2/(rr+1)^2]))) + log((rr+1)/(rr-1))/rr*ellipke((rr^2-1)/rr^2) );
      C_mm_h_DC = epsilon_mm / s * t*Length_out(j);
%      C_mm_h_DC_f = C_mm_h_DC - epsilon_mm_eq / s * t*Length_out(j);
      if (j == i+1)
        C_mm_h = C_mm_h +  C_mm_h_DC * (sum_r(j)-sum_r(i)+(r(i)-r(j))/2)^2;  %portion of each turn
      else
%        C_mm_h = C_mm_h +  0.5*C_mm_h_DC_f* (sum_r(j)-sum_r(i)+(r(i)-r(j))/2)^2;  %portion of each turn
      end
      %C_mm_h = C_mm_h +  (C_mm_h_DC_f-C_mm_h_DC_a/epsilon_mm*epsilon_mmf+C_mm_h_DC_a)* (sum_r(j)-sum_r(i)+(r(i)-r(j))/2)^2;  %portion of each turn
      %C_mm_h = C_mm_h + C_mm_h_DC_f*(sum_r(j)-sum_r(i)+(r(i)-r(j))/2)^2;
    end
  end
  C_substrate;
  C_mm_h;
  C_mm_v    ;
  % Total effective capacitance is sum of the capacitance between metal and substrate and the capacitance between metals of adjacent turns
  C = C_substrate + C_mm_h + C_mm_v;
  
end