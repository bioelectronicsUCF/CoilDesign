function [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t)
  
   u = pi*(4e-7)*1e-6;     % Permeability of free space.
   mu = u;


   I_eddy_edge = linspace(0,0,n);  
   I_eddy_center = linspace(0,0,n);
   M_eddy_ec = u * pi * (3/8*w)^2 / (2*((w/2-delta/2)^2+(3/8*w)^2)^(3/2))*(1+15/32*(2*3/8*w*(w/2-delta/2)/((3/8*w)^2+(w/2-delta/2)^2))^2+315/1024*(2*3/8*w*(w/2-delta/2)/((3/8*w)^2+(w/2-delta/2)^2))^4);sM_=1i*omega*M_eddy_ec;
   V_eddy_edge = double(1i * omega * mu * (H * (w-2*delta))*1e6);              V_e=double(V_eddy_edge);
   V_eddy_center = double(1i * omega * mu * H * w/2*1e6);                      V_c=V_eddy_center;
   if (w > 5*delta)
     if (t > 2*delta)
       R_eddy_edge = 2*Sheet_R_coil/delta;                         R_e=R_eddy_edge;
       R_eddy_center = 2*(Sheet_R_coil*t/(2*delta))/(w/4-delta);   R_c=R_eddy_center;
       L_eddy_edge = mu/pi*log((w-delta)/delta);                   sL_e=1i*omega*L_eddy_edge;
       L_eddy_center = mu/pi*log((0.75*w-delta)/(0.25*w-delta));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = (V_e-sM_/(R_c+sL_c)*V_c) / (R_e+sL_e-sM_^2/(R_c+sL_c)) /1e6;
       I_eddy_center = (V_c-sM_/(R_e+sL_e)*V_e) / (R_c+sL_c-sM_^2/(R_e+sL_e))/1e6;
     else
       R_eddy_center = 2*Sheet_R_coil/(w/4);   R_c=R_eddy_center;
       L_eddy_center = mu/pi*log(w/(w/4));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = linspace(0,0,n);
       I_eddy_center = V_c / (R_c+sL_c)/1e6;
     end
   elseif (w > 2*delta)
     R_eddy_edge = 2*Sheet_R_coil/delta;                         R_e=R_eddy_edge;
     L_eddy_edge = mu/pi*log((w-delta)/delta);                   sL_e=1i*omega*L_eddy_edge;
     I_eddy_edge = V_e / (R_e+sL_e)/1e6;
     I_eddy_center = linspace(0,0,n);
   else
     if (t > 2*delta)
       R_eddy_center = 2*(Sheet_R_coil*t/(2*delta))/(w/4-delta);   R_c=R_eddy_center;
       L_eddy_center = mu/pi*log((0.75*w-delta)/(0.25*w-delta));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = linspace(0,0,n);
       I_eddy_center = V_c / (R_c+sL_c)/1e6;
     else
       R_eddy_center = 2*Sheet_R_coil/(w/4);   R_c=R_eddy_center;
       L_eddy_center = mu/pi*log(0.75*w/(0.25*w));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = linspace(0,0,n);
       I_eddy_center = V_c / (R_c+sL_c)/1e6;
     end
   end

end