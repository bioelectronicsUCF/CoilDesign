
function [I_eddy_edge, I_eddy_center] = Ieddy7(w,delta,H,Sheet_R_coil,omega,n,t,length_turn)
  
   u = pi*(4e-7)*1e-6;     % Permeability of free space.
   mu = u;

   I_eddy_edge = linspace(0,0,n);  
   I_eddy_center = linspace(0,0,n);
   
%   M_eddy_ec = u * pi * length_turn * (3/8*w)^2 / (2*((w/2-delta/2)^2+(3/8*w)^2)^(3/2))*(1+15/32*(2*3/8*w*(w/2-delta/2)/((3/8*w)^2+(w/2-delta/2)^2))^2+315/1024*(2*3/8*w*(w/2-delta/2)/((3/8*w)^2+(w/2-delta/2)^2))^4);sM_=1i*omega*M_eddy_ec;
   M_eddy_ec = mu/pi.*length_turn.*(log((1+length_turn.^2/w^2).^0.5+length_turn/w)-(1+w^2./length_turn.^2).^0.5+w./length_turn)-mu/pi*length_turn.*(log((1+(0.25*length_turn).^2/(0.75*w-delta)^2).^0.5+(0.25*length_turn)/(0.75*w-delta))-(1+(0.75*w-delta)^2./(0.25*length_turn).^2).^0.5+(0.75*w-delta)./(0.25*length_turn));
   sM_=1i*omega*M_eddy_ec;
   V_eddy_edge = double(1i * omega * mu * (H .* length_turn * (w-2*delta))*1e6);              V_e=double(V_eddy_edge);
   V_eddy_center = double(1i * omega * mu * H .* length_turn * w/2*1e6);                      V_c=V_eddy_center;
   if (w > 2*delta)
     if (t > 2*delta)
       R_eddy_edge = 2*Sheet_R_coil * length_turn /delta;                         R_e=R_eddy_edge;
       R_eddy_center = 2*(Sheet_R_coil*t/(2*delta)) * length_turn /(w/2-delta);   R_c=R_eddy_center;
       L_eddy_edge = mu/pi * length_turn *log(4*(w-delta)/(delta+t));                   sL_e=1i*omega*L_eddy_edge;
       L_eddy_center = mu/pi * length_turn *log((2*w-4*delta)/(0.5*w+delta));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = (V_e-sM_./(R_c+sL_c).*V_c) ./ (R_e+sL_e-sM_.^2./(R_c+sL_c)) /1e6;
       I_eddy_center = (V_c-sM_./(R_e+sL_e).*V_e) ./ (R_c+sL_c-sM_.^2./(R_e+sL_e))/1e6;
     else
       R_eddy_edge = 2*Sheet_R_coil * length_turn /delta;                         R_e=R_eddy_edge;
       R_eddy_center = 2*Sheet_R_coil * t/delta * length_turn  /(1-exp(-t/delta))/(w/2);   R_c=R_eddy_center;
       L_eddy_edge = mu/pi * length_turn *log((w-delta)/delta);                   sL_e=1i*omega*L_eddy_edge;
       L_eddy_center = mu/pi * length_turn *log(w/(w/2));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = (V_e-sM_./(R_c+sL_c).*V_c) ./ (R_e+sL_e-sM_.^2./(R_c+sL_c)) /1e6;
       I_eddy_center = (V_c-sM_./(R_e+sL_e).*V_e) ./ (R_c+sL_c-sM_.^2./(R_e+sL_e))/1e6;
     end
   else
     if (t > 2*delta)
       R_eddy_center = 2*(Sheet_R_coil*t/(2*delta)) * length_turn /(w/2);   R_c=R_eddy_center;
       L_eddy_center = mu/pi * length_turn *log((4*w-8*delta)/(0.5*w+delta));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = linspace(0,0,n);
       I_eddy_center = V_c ./ (R_c+sL_c)/1e6;
     else
       R_eddy_center = 2*Sheet_R_coil * length_turn /(w/4);   R_c=R_eddy_center;
       L_eddy_center = mu/pi * length_turn *log(0.75*w/(0.25*w));     sL_c=1i*omega*L_eddy_center;
       I_eddy_edge = linspace(0,0,n);
       I_eddy_center = V_c ./ (R_c+sL_c)/1e6;
     end
   end

end

