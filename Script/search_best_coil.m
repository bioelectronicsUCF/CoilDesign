% Renamed from the original: "SquareSpiralOptimization.m" in the Lake Nona Matlab folder.
% Input parameters: distance between the coils (d), outer diameter (OD) and
% inner diameter (ID) of the coils, manufacturing technology data, operating frequency, 
% target load impedance, and expected access resistance.
% Presets are provided for manufacturing technology data. 
% User can choose to adjust each manufacturing technology data to improve the accuracy of PTE estimation.
% Output parameters: maximum achievable PTE, specifications of the best coil set. 
clear all 
close all
clc

%% Input parameters:
%Distance between the coils (d), outer diameter (OD) and inner diameter (ID) of the coils
d = 15000.0;                                %in um, %Distance between the Transmitting and Receiving coils
Dout_RX = 15000.0;                          %in um, %Outer diameter of the Receiving Coil 
Din_RX = 10000.0;                           %in um, %Inner diameter of the Receiving Coil
Dout_TX = 50000.0;                          %in um, %Outer diameter of the Transmitting Coil 
Din_TX = 20000.0;                           %in um, %Inner diameter of the Transmitting Coil
%Preset for manufacturing technology data of RX coil: MEMS process on custom IC 
max_thickness_RX = 30.0;                    %in um, %Maximum thickness of the conductor traces
min_width_RX = 10.0;                        %in um, %Minimum width of conducting traces
min_space_RX = 10.0;                        %in um, %Minimum spacing between conducting traces
rho_RX = 2.65e-2;                           %in ohm*um, %Resistance of the conducting material
epsilon_RX = 3.63e-17;                      %in Farad/um, %Permittivity of the insulating material
%Preset for manufacturing technology data of TX coil: PCB process
max_thickness_TX = 500.0;                   %in um, %Maximum thickness of the conductor traces
min_width_TX = 100.0;                       %in um, %Minimum width of conducting traces
min_space_TX = 100.0;                       %in um, %Minimum spacing between conducting traces
rho_TX = 1.72e-2;                           %in ohm*um, %Resistance of the conducting material
epsilon_TX = 3.10e-17 - 1j * 2.48e-19;      %in Farad/um, %Permittivity of the insulating material
%Application-specific constraints
f = 13.56e6;                                %in Hz, %Operating frequency of the wireless power transfer
RL = 100;                                   %in ohm, %Target load resistance
RA = 0.20;                                  %in ohm, %Expected access resistance
%Optimization-related constraints
n_search_points = 10;                       %number of points corresponding each variable per each iteration to calculate PTE 
iteration = 1;                              %number of iterations for optimization
max_turn_coil = 15;                         %maximum number of turns of coils to exploit

%% Generating a dummy set of TX and RX coils:
thickness_RX = max_thickness_RX;
turn_RX = 3;
space_RX = min_space_RX;
width_RX = ((Dout_RX-Din_RX)/2+space_RX)/turn_RX-space_RX;
sheet_RX = rho_RX / thickness_RX;

R_RX = Res_tot6(space_RX, turn_RX, Dout_RX, Din_RX, sheet_RX, sheet_RX, 0, thickness_RX, f);
L_RX = wheelerTurns(space_RX, turn_RX, Dout_RX, Din_RX)   ;
C_RX = distributedCap5(space_RX, turn_RX, Dout_RX, Din_RX, epsilon_RX, epsilon_RX, epsilon_RX, 5, 10, thickness_RX, 0, 0, f);

thickness_TX = max_thickness_TX;
turn_TX = 5;
space_TX = (Dout_TX-Din_TX)/(2*turn_TX+2*(turn_TX-1));
width_TX = ((Dout_TX-Din_TX)/2+space_TX)/turn_TX-space_TX;
sheet_TX = rho_TX / thickness_TX;

R_TX = RA+Res_tot6(space_TX, turn_TX, Dout_TX, Din_TX, sheet_TX, sheet_TX, 0, thickness_TX, f);
L_TX = wheelerTurns(space_TX, turn_TX, Dout_TX, Din_TX)   ;
C_TX = distributedCap5(space_TX, turn_TX, Dout_TX, Din_TX, epsilon_TX, epsilon_TX, epsilon_TX, 1000, 1000, thickness_TX, 0, 0, f);

M = MFilament2(space_RX, turn_RX, Dout_RX, Din_RX, space_TX, turn_TX, Dout_TX, Din_TX, d+thickness_TX/2+thickness_RX/2);

%% Start iteration on RX coils
%check input parameters
max_turn_RX_possible=((Dout_RX-Din_RX)/2+min_space_RX)/(min_width_RX+min_space_RX);
if max_turn_coil > max_turn_RX_possible
    max_turn_RX = max_turn_RX_possible;
else
    max_turn_RX = max_turn_coil;
end

%setup iteration range
thickness_RX_vals = linspace(max_thickness_RX/n_search_points,max_thickness_RX,n_search_points);
turn_RX_vals = round(linspace(1,max_turn_RX,n_search_points));
space_RX_vals = zeros(n_search_points,n_search_points);
max_space_RX = zeros(n_search_points,1);
eta_RX = zeros(n_search_points,n_search_points,n_search_points);
for i_turn = 1 : n_search_points
    turn_RX = turn_RX_vals(i_turn);
    if turn_RX == 1
        space_RX_vals(i_turn,1:n_search_points) = linspace(Din_RX,Din_RX,n_search_points);
    else
        max_space_RX(i_turn) = ((Dout_RX-Din_RX)/2-min_width_RX)/(turn_RX-1)-min_width_RX;
        space_RX_vals(i_turn,1:n_search_points) = linspace(min_space_RX,max_space_RX(i_turn),n_search_points);
    end
end

max_eta=0;
%calculate PTE for all virtual RX coils
for i_thickness = 1:n_search_points
    for i_turn = 1:n_search_points
        for i_space = 1:n_search_points
            thickness_RX = thickness_RX_vals(i_thickness);
            turn_RX = turn_RX_vals(i_turn);
            space_RX = space_RX_vals(i_turn,i_space);
            width_RX = ((Dout_RX-Din_RX)/2+space_RX)/turn_RX-space_RX;
            sheet_RX = rho_RX / thickness_RX;
            R_RX = Res_tot6(space_RX, turn_RX, Dout_RX, Din_RX, sheet_RX, sheet_RX, 0, thickness_RX, f);
            L_RX = wheelerTurns(space_RX, turn_RX, Dout_RX, Din_RX)   ;
            C_RX = distributedCap5(space_RX, turn_RX, Dout_RX, Din_RX, epsilon_RX, epsilon_RX, epsilon_RX, 5, 10, thickness_RX, 0, 0, f);
            M = MFilament2(space_RX, turn_RX, Dout_RX, Din_RX, space_TX, turn_TX, Dout_TX, Din_TX, d+thickness_TX/2+thickness_RX/2);
            eta_RX(i_thickness,i_turn,i_space) = maxPTEc_sweep_cap3(R_TX,L_TX,C_TX,R_RX,L_RX,C_RX,M,f,RL,5);
            if eta_RX(i_thickness,i_turn,i_space) > max_eta
                max_eta = eta_RX(i_thickness,i_turn,i_space);
                output_thickness_RX = thickness_RX;
                output_turn_RX = turn_RX;
                output_space_RX = space_RX;
            end
        end
    end
end

%Show optimum thickness, turn count, spacing, and width of RX coil on
%command line
output_thickness_RX
output_turn_RX
output_space_RX
output_width_RX = ((Dout_RX-Din_RX)/2+output_space_RX)/output_turn_RX-output_space_RX
%

%% Start iteration on TX coils
%Update RX coil parameters with the optimum values obtained above
thickness_RX = output_thickness_RX;
turn_RX = output_turn_RX;
space_RX = output_space_RX;
width_RX = ((Dout_RX-Din_RX)/2+space_RX)/turn_RX-space_RX;
sheet_RX = rho_RX / thickness_RX;

R_RX = Res_tot6(space_RX, turn_RX, Dout_RX, Din_RX, sheet_RX, sheet_RX, 0, thickness_RX, f);
L_RX = wheelerTurns(space_RX, turn_RX, Dout_RX, Din_RX)   ;
C_RX = distributedCap5(space_RX, turn_RX, Dout_RX, Din_RX, epsilon_RX, epsilon_RX, epsilon_RX, 5, 10, thickness_RX, 0, 0, f);


%check input parameters
max_turn_TX_possible=((Dout_TX-Din_TX)/2+min_space_TX)/(min_width_TX+min_space_TX);
if max_turn_coil > max_turn_TX_possible
    max_turn_TX = max_turn_TX_possible;
else
    max_turn_TX = max_turn_coil;
end

%setup iteration range
thickness_TX_vals = linspace(max_thickness_TX/n_search_points,max_thickness_TX,n_search_points);
turn_TX_vals = round(linspace(1,max_turn_TX,n_search_points));
space_TX_vals = zeros(n_search_points,n_search_points);
max_space_TX = zeros(n_search_points,1);
eta_TX = zeros(n_search_points,n_search_points,n_search_points);
for i_turn = 1 : n_search_points
    turn_TX = turn_TX_vals(i_turn);
    if turn_TX == 1
        space_TX_vals(i_turn,1:n_search_points) = linspace(Din_TX,Din_TX,n_search_points);
    else
        max_space_TX(i_turn) = ((Dout_TX-Din_TX)/2-min_width_TX)/(turn_TX-1)-min_width_TX;
        space_TX_vals(i_turn,1:n_search_points) = linspace(min_space_TX,max_space_TX(i_turn),n_search_points);
    end
end

max_eta = 0;
for i_thickness = 1:n_search_points
    for i_turn = 1:n_search_points
        for i_space = 1:n_search_points
            thickness_TX = thickness_TX_vals(i_thickness);
            turn_TX = turn_TX_vals(i_turn);
            space_TX = space_TX_vals(i_turn,i_space);
            width_TX = ((Dout_TX-Din_TX)/2+space_TX)/turn_TX-space_TX;
            sheet_TX = rho_TX / thickness_TX;
            R_TX = RA+Res_tot6(space_TX, turn_TX, Dout_TX, Din_TX, sheet_TX, sheet_TX, 0, thickness_TX, f);
            L_TX = wheelerTurns(space_TX, turn_TX, Dout_TX, Din_TX)   ;
            C_TX = distributedCap5(space_TX, turn_TX, Dout_TX, Din_TX, epsilon_TX, epsilon_TX, epsilon_TX, 5, 10, thickness_TX, 0, 0, f);
            M = MFilament2(space_RX, turn_RX, Dout_RX, Din_RX, space_TX, turn_TX, Dout_TX, Din_TX, d+thickness_TX/2+thickness_RX/2);
            eta_TX(i_thickness,i_turn,i_space) = maxPTEc_sweep_cap3(R_TX,L_TX,C_TX,R_RX,L_RX,C_RX,M,f,RL,5);
            if eta_TX(i_thickness,i_turn,i_space) > max_eta
                max_eta = eta_TX(i_thickness,i_turn,i_space);
                output_thickness_TX = thickness_TX;
                output_turn_TX = turn_TX;
                output_space_TX = space_TX;
            end
        end
    end
end

%Show optimum thickness, turn count, spacing, and width of TX coil on
%command line
output_thickness_TX
output_turn_TX
output_space_TX
output_width_TX = ((Dout_TX-Din_TX)/2+output_space_TX)/output_turn_TX-output_space_TX
%

%% Visualize PTE with respect to parameters
%Create window to draw figures
figure('position',[100 0 300 450])


%Plot PTE vs TX parameters
subplot(3,2,2)
pte_plot3(thickness_TX_vals,100*max(eta_TX,[],[2,3]),"thickness_{TX} ({\mu}m)","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))
[max_eta,i_output_thickness_TX]=max(max(eta_TX,[],[2,3]));

subplot(3,2,4)
pte_plot3(turn_TX_vals,reshape(100*max(eta_TX(i_output_thickness_TX,:,:),[],[3]),1,[]),"turn count_{TX}","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))
[max_eta,i_output_turn_TX]=max(max(eta_TX(i_output_thickness_TX,:,:),[],[3]));

subplot(3,2,6)
pte_plot3(reshape(space_TX_vals(i_output_turn_TX,:),1,[]),reshape(100*eta_TX(i_output_thickness_TX,i_output_turn_TX,:),1,[]),"spacing_{TX} {\mu}m","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))


% Generating the best TX and RX coils:
thickness_RX = output_thickness_RX;
turn_RX = output_turn_RX;
space_RX = output_space_RX;
width_RX = ((Dout_RX-Din_RX)/2+output_space_RX)/output_turn_RX-output_space_RX;
sheet_RX = rho_RX / output_thickness_RX;

R_RX = Res_tot6(space_RX, turn_RX, Dout_RX, Din_RX, sheet_RX, sheet_RX, 0, thickness_RX, f);
L_RX = wheelerTurns(space_RX, turn_RX, Dout_RX, Din_RX)   ;
C_RX = distributedCap5(space_RX, turn_RX, Dout_RX, Din_RX, epsilon_RX, epsilon_RX, epsilon_RX, 5, 10, thickness_RX, 0, 0, f);

thickness_TX = output_thickness_TX;
turn_TX = 5;
space_TX = (Dout_TX-Din_TX)/(2*turn_TX+2*(turn_TX-1));
width_TX = ((Dout_TX-Din_TX)/2+space_TX)/turn_TX-space_TX;
sheet_TX = rho_TX / thickness_TX;

R_TX = RA+Res_tot6(space_TX, turn_TX, Dout_TX, Din_TX, sheet_TX, sheet_TX, 0, thickness_TX, f);
L_TX = wheelerTurns(space_TX, turn_TX, Dout_TX, Din_TX)   ;
C_TX = distributedCap5(space_TX, turn_TX, Dout_TX, Din_TX, epsilon_TX, epsilon_TX, epsilon_TX, 1000, 1000, thickness_TX, 0, 0, f);

M = MFilament2(space_RX, turn_RX, Dout_RX, Din_RX, space_TX, turn_TX, Dout_TX, Din_TX, d+thickness_TX/2+thickness_RX/2);

max_eta = 0;
%Calculate and plot PTE vs RX parameters
for i_thickness = 1:n_search_points
    for i_turn = 1:n_search_points
        for i_space = 1:n_search_points
            thickness_RX = thickness_RX_vals(i_thickness);
            turn_RX = turn_RX_vals(i_turn);
            space_RX = space_RX_vals(i_turn,i_space);
            width_RX = ((Dout_RX-Din_RX)/2+space_RX)/turn_RX-space_RX;
            sheet_RX = rho_RX / thickness_RX;
            R_RX = Res_tot6(space_RX, turn_RX, Dout_RX, Din_RX, sheet_RX, sheet_RX, 0, thickness_RX, f);
            L_RX = wheelerTurns(space_RX, turn_RX, Dout_RX, Din_RX)   ;
            C_RX = distributedCap5(space_RX, turn_RX, Dout_RX, Din_RX, epsilon_RX, epsilon_RX, epsilon_RX, 5, 10, thickness_RX, 0, 0, f);
            M = MFilament2(space_RX, turn_RX, Dout_RX, Din_RX, space_TX, turn_TX, Dout_TX, Din_TX, d+thickness_TX/2+thickness_RX/2);
            eta_RX(i_thickness,i_turn,i_space) = maxPTEc_sweep_cap3(R_TX,L_TX,C_TX,R_RX,L_RX,C_RX,M,f,RL,5);
            if eta_RX(i_thickness,i_turn,i_space) > max_eta
                max_eta = eta_RX(i_thickness,i_turn,i_space);
                output_thickness_RX = thickness_RX;
                output_turn_RX = turn_RX;
                output_space_RX = space_RX;
            end
        end
    end
end

%Plot PTE vs RX parameters
subplot(3,2,1)
pte_plot3(thickness_RX_vals,100*max(eta_RX,[],[2,3]),"thickness_{RX} ({\mu}m)","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))
[max_eta,i_output_thickness_RX]=max(max(eta_RX,[],[2,3]));

subplot(3,2,3)
pte_plot3(turn_RX_vals,reshape(100*max(eta_RX(i_output_thickness_RX,:,:),[],[3]),1,[]),"turn count_{RX}","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))
[max_eta,i_output_turn_RX]=max(max(eta_RX(i_output_thickness_RX,:,:),[],[3]));

subplot(3,2,5)
pte_plot3(reshape(space_RX_vals(i_output_turn_RX,:),1,[]),reshape(100*eta_RX(i_output_thickness_RX,i_output_turn_RX,:),1,[]),"spacing_{RX} {\mu}m","PTE (%)",0,'a','a',0,10*ceil(max_eta*10),ceil(max_eta*10))


%% Define functions:
% Given a spacing, coil diameters, and integer of complete turns; return
% Exp[L] and arbitrarily precise trace Width.
function [L, w] = wheelerTurns(s, n, Dout, Din)
    
   % Useful definitions to carry through:
   delta = (Dout - Din)/2;
   Davg  = (Dout + Din)/2;
   rho   = (Dout - Din) / (Dout + Din);
   
   % Calculate the width from the given (integer) number of turns:
   w = (delta - s*(n-1))/n;
   
   % Using Modified Wheeler approximation from Mohan's 1999 paper.
   K1 = 2.34;       % Table 1: K1 Square.
   K2 = 2.75;       % Table 1: K2 Square.
   u = pi*(4e-7);     % Permeability of free space.
   
   L = K1*u*n*n*Davg/(1+K2*rho)*1e-6;

end


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


function Rprox_tot = Res_prox_only15(s, n, Dout, Din, Sheet_R_coil, Sheet_R_con, l_con, t, f)

  w = ((Dout - Din)/2 - s*(n-1))/n;
  length = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
  
  rho = Sheet_R_coil*t;
  mu = 4 * pi * 1e-7 * 1e-6;
  delta = (rho/pi/f/mu)^0.5;
  
  omega = 2 * pi * f;
  
  for i = 1:n
    length_turn(i) = (4*(Din+w)+8*(i-1)*(w+s));
    R_dc(i) = Sheet_R_coil / w * length_turn(i);
  end

  
  k = (1+1j)/delta;
  
  Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );

  Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
  Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
  Hin_ = Hin + linspace(0,(n-1)/n,n) * (Hin-Hout);
  Hout_ = Hout - linspace((n-1)/n,0,n) * (Hout-Hin);
  H = Hin + (linspace(0,(n-1)/n,n) + 1/2/n) * (Hin-Hout);
  
  delH = (Hin-Hout)/n;
  delH_ = linspace(1,1,n) * delH;
  
  
  for i = 1:n
    len_ref = Din + w + 2*(i-1)*(w+s);

    Hin_(i)=0;

    b2 = len_ref/2;
    b3 = len_ref-0.5*w;
    b4 = len_ref/2;  
    Hin_(i)=Hin_(i)+1/4/pi/b2*b3/sqrt(b2^2+b3^2)+1/4/pi/b3*(b2/sqrt(b2^2+b3^2)+b4/sqrt(b3^2+b4^2))+1/4/pi/b4*b3/(b3^2+b4^2);
%    Hin_(i)=Hin_(i)+1/2/pi/w*abs(1+0.5*log(1+(w/t)^2)-w/t*atan(t/w));

    for j = i+1 : n            
      c1 = 0.5*w+(j-i)*(w+s);
      c2 = len_ref/2+(j-i)*(w+s);
      c3 = len_ref-0.5*w+(j-i)*(w+s);
      c4 = len_ref/2+(j-i)*(w+s);
      Hin_(i)=Hin_(i)+1/4/pi/c1*(c4/sqrt(c1^2+c4^2)+c2/sqrt(c1^2+c2^2))+1/4/pi/c2*(c1/sqrt(c1^2+c2^2)+c3/sqrt(c2^2+c3^2))+1/4/pi/c3*(c2/sqrt(c2^2+c3^2)+c4/sqrt(c3^2+c4^2))+1/4/pi/c4*(c3/sqrt(c3^2+c4^2)+c1/sqrt(c1^2+c4^2));
    end

    for j = 1 : i-1            
      a1 = -0.5*w+(i-j)*(w+s);
      a2 = len_ref/2-(i-j)*(w+s);
      a3 = len_ref-0.5*w-(i-j)*(w+s);
      a4 = len_ref/2-(i-j)*(w+s);
      Hin_(i)=Hin_(i)-1/4/pi/a1*(a4/sqrt(a1^2+a4^2)+a2/sqrt(a1^2+a2^2))+1/4/pi/a2*(a1/sqrt(a1^2+a2^2)+a3/sqrt(a2^2+a3^2))+1/4/pi/a3*(a2/sqrt(a2^2+a3^2)+a4/sqrt(a3^2+a4^2))+1/4/pi/a4*(a3/sqrt(a3^2+a4^2)+a1/sqrt(a1^2+a4^2));
    end
    
    Hout_(i)=0;
    for j = 1 : i-1            
      a1 = 0.5*w+(i-j)*(w+s);
      a2 = len_ref/2-(i-j)*(w+s);
      a3 = len_ref+0.5*w-(i-j)*(w+s);
      a4 = len_ref/2-(i-j)*(w+s);
      Hout_(i)=Hout_(i)-1/4/pi/a1*(a4/sqrt(a1^2+a4^2)+a2/sqrt(a1^2+a2^2))+1/4/pi/a2*(a1/sqrt(a1^2+a2^2)+a3/sqrt(a2^2+a3^2))+1/4/pi/a3*(a2/sqrt(a2^2+a3^2)+a4/sqrt(a3^2+a4^2))+1/4/pi/a4*(a3/sqrt(a3^2+a4^2)+a1/sqrt(a1^2+a4^2));
      -1/4/pi/a1*(a4/sqrt(a1^2+a4^2));
    end

    for j = i+1 : n            
      c1 = -0.5*w+(j-i)*(w+s);
      c2 = len_ref/2+(j-i)*(w+s);
      c3 = len_ref+0.5*w+(j-i)*(w+s);
      c4 = len_ref/2+(j-i)*(w+s);
      Hout_(i)=Hout_(i)+1/4/pi/c1*(c4/sqrt(c1^2+c4^2)+c2/sqrt(c1^2+c2^2))+1/4/pi/c2*(c1/sqrt(c1^2+c2^2)+c3/sqrt(c2^2+c3^2))+1/4/pi/c3*(c2/sqrt(c2^2+c3^2)+c4/sqrt(c3^2+c4^2))+1/4/pi/c4*(c3/sqrt(c3^2+c4^2)+c1/sqrt(c1^2+c4^2));
    end

    b2 = len_ref/2;
    b3 = len_ref-0.5*w;
    b4 = len_ref/2;  
    Hout_(i)=Hout_(i)+1/4/pi/b2*b3/sqrt(b2^2+b3^2)+1/4/pi/b3*(b2/sqrt(b2^2+b3^2)+b4/sqrt(b3^2+b4^2))+1/4/pi/b4*b3/(b3^2+b4^2);
%    Hout_(i)=Hout_(i)-1/2/pi/w*abs(1+0.5*log(1+(w/t)^2)-w/t*atan(t/w));
  end
  H=(Hin_+Hout_)/2;
  delH_=Hin_-Hout_  ;

   

  [I_eddy_edge, I_eddy_center] = Ieddy7(w,delta,H,Sheet_R_coil,omega,n,t,length_turn);
  
  delH_inc = linspace(0,0,n);
  Hin_inc = linspace(0,0,n);
  Hout_inc = linspace(0,0,n);
  
%{  
  for i = 1:n
    for j = 1:n
      if (j==i)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(t+delta*2)-1/(w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(s+delta/2)-2/(s+w-delta/2)+1/(s+2*w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(w/8+s)-1/(s+w*7/8)-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
%        delH_inc(i) = delH_inc(i) + I_eddy_edge(j)/2/pi * (1/(s+delta/2)-2/(s+w-delta/2)+1/(s+2*w-delta/2));
%        delH_inc(i) = delH_inc(i) + I_eddy_center(j)/2/pi * (1/(w/8+s)-1/(s+w*7/8)-1/(s+w*9/8)+1/(s+15/8*w));
      end
    end
  end
  
  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/5/pi * (1/(t+delta*2)-1/(w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/5/pi * (1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/4/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      elseif (j==i+1)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)/4/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)/4/pi * (1/(w/8+s)-1/(s+w*7/8));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)/4/pi * (-1/(s+w-delta/2)+1/(s+2*w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)/4/pi * (-1/(s+w*9/8)+1/(s+15/8*w));
      end
    end
  end
%}

  for i = 1:n
    for j = 1:n
      if (abs(i-j)==0)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j) * (1/(t*2+delta*2)-1/(w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j) * 2/pi *(1/(t+w/2)-1/(w-w/4));
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j) * (1/(t*2+delta*2)-1/(w-delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j) * 2/pi *(1/(t+w/2)-1/(w-w/4));
      elseif (j==i-1)
        Hout_inc(i) = Hout_inc(i) - I_eddy_edge(j)*2/pi * (1/(s+delta/2));
        Hout_inc(i) = Hout_inc(i) - I_eddy_center(j)*2/pi * (1/(w/8+s));
      elseif (j==i+1)
        Hin_inc(i) = Hin_inc(i) + I_eddy_edge(j)*2/pi * (1/(s+delta/2)-1/(s+w-delta/2));
        Hin_inc(i) = Hin_inc(i) + I_eddy_center(j)*2/pi * (1/(w/8+s)-1/(s+w*7/8));
      end
    end
  end

%  delH_ = delH_ + delH_inc;
  Hin_ = Hin_+Hin_inc;
  Hout_ = Hout_+Hout_inc;
  delH_ = Hin_ - Hout_;


  HE_ = (1-delH_*t)/2/w;  
%  HE = 1*(1-delH*t)/2/w;
%  HE_ = linspace(1,1,n) * HE;
  
  HE_inc = linspace(0,0,n);

  HE_ = HE_ + HE_inc;
  
  alpha = HE_./delH_;
  beta = Hin./delH_;
  
  Rprox_tot=0;
  Rdc_tot=0;
  for i = 1:n
    cothkw=coth(k*w);
    if(isnan(cothkw))
      if(real(k*w)>0)
        cothkw=1;
      else
        cothkw=-1;
      end
    end
    cothtk=coth(t*k/2);
    if(isnan(cothtk))
      if(real(t*k/2)>0)
        cothtk=1;
      else
        cothtk=-1;
      end
    end
    tanhkw=tanh(k*w/2);
    if(isnan(tanhkw))
      if(real(k*w)>0)
        tanhkw=1;
      else
        tanhkw=-1;
      end
    end
    Rprox_per_dc(i) = abs(real(w*t/(2*alpha(i)*w+t)^2*(t*k*cothkw+2*k*w*alpha(i)^2*cothtk+4*alpha(i)+2*t*k*tanhkw*abs(Hin_(i))*abs(Hout_(i))/abs(delH_(i))^2)));%*(beta(i)+1-i)*(beta(i)-i))));%
    Rprox(i) = R_dc(i)*Rprox_per_dc(i);
    Rprox_tot = Rprox_tot + Rprox(i);
    Rdc_tot = Rdc_tot + R_dc(i);
  end
end

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

% Given a spacing, coil diameters, turn count, insulator relative permittivity,
% insulator height, metal thickness; return effective capacitance
function [C, C_substrate, C_mm_h, C_mm_v] = distributedCap5(s, n, od, id, epsilon_r_sub, epsilon_r_mm, epsilon_r_mmf, d, h, t, c_sub_per_area, g_sub_per_area, f)


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
    if n > 15
      for i = 1 : n
        C_turns_sub(i) = epsilon_sub / h * Area(i);
        C_sub(i) = Area(i) * c_sub_per_area;
        R_sub(i) = 1 / Area(i) / g_sub_per_area;
        C_substrate = C_substrate + b(i) * 1/( 1/C_turns_sub(i)+(2i*pi*f*R_sub(1))/(2i*pi*f*R_sub(1)*C_sub(i)+1) );
      end
    else
      Area_tot = sum(Area);
      C_turns_sub_tot = epsilon_sub / h * Area_tot;
      C_sub_tot = c_sub_per_area * Area_tot;
      g_sub_tot = g_sub_per_area * Area_tot;
      C_substrate = 1/6 * 1/(1/C_turns_sub_tot+1/(C_sub_tot+g_sub_tot/2i/pi/f));
    end
  end
    
  if n > 0
    C_mm_v = epsilon_mm / d * w^2 * b(1);
    for i = 2 : n
      C_mm_v = C_mm_v + epsilon_mm / d * w^2 * b(i);
    end
  else
    C_mm_v = 0;
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
%     rr = (s*abs(j-i)+w*2/3)/(s*abs(j-i));
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

function M = MFilament2(rxsp, rxn, rxod, rxid, txsp, txn, txod, txid, z)

  txw = ((txod-txid)/2+txsp)/txn-txsp;
  rxw = ((rxod-rxid)/2+rxsp)/rxn-rxsp;

  rho = (4/pi)^(1+rxod/txod);
  mu = 4 * pi * 1e-7 * 1e-6;  %H/um
  
  M = 0;
  for tx_seg = 1 : txn
    for rx_seg = 1 : rxn
      a = txod/2 - (tx_seg-1)*(txw+txsp) - txw/2;
      b = rxod/2 - (rx_seg-1)*(rxw+rxsp) - rxw/2;
      gamma = 2 * a * b / (a^2 + b^2 + z^2);
      
      M_ = mu * pi * a^2 * b^2 / (2 * (a^2 + b^2 + z^2) ^ 1.5) * (1 + 15/32 * gamma^2 + 315/1024 *gamma^4);
      M = M + rho * M_;
    end
  end

end
% Input: effective R and L of the TX and RX coils, tuning cap sweep range and resolution, k between TX and RX coils, frequency, load resistance, and target received voltage.
% Output: maximum power transfer efficiency, tuning caps, voltages and currents of TX and RX ports
%
%
% Circuit model
%
%     VT                 VM                         VR
%     *-||-----^^^----mmm-*-mmm---^^^---------------*
%       CT  |  RT    LT-M | LR-M   RR  |      |     |
%          ---            3           ---    ---    <
%      CPT ---          M 3       CPR --- CL --- RL <   
%           |             |            |      |     |
%     -----------------------------------------------
%
% 
% The circuit will be simplified as below before calculating power transfer efficiency
%
%     VT       VTin         VM             VR
%     *----||-------|ZT|----*-----|ZR|-----*
%     -->  CT   |   -->     |              |       
%     IT       ---  ITin    |              |
%          CPT ---        |sM|           |ZL|
%               |           |              | 
%     --------------------------------------
%

function [eta,CT,CL,VT,VM,IT] = maxPTEc_sweep_cap3(RT,LT,CPT,RR,LR,CPR,M,f,RL,VR)

s = 1i * 2 * pi * f; sLT = s * LT; sLR = s * LR; sM = s * M;      %Impedances of inductances at frequency f    
w = 2 * pi * f;
          
%CT = real(-1/s / (1/(s*CPT+1/(RT+sLT))));

CL_max = 1/( RR^2/LR + w^2*LR );
CL_init = real(CPR);
if CL_max > CL_init
  CL = CL_max - CL_init;
else
  CL = 0;
end

  
  sCR = s*(CL+CPR);
  
  ZR = RR + sLR - sM;                                       %ZR is calculated by adding impedances of RR and LR
  ZL = 1 / ( sCR + 1 / RL );                                %ZL is calculated by adding admittances of CPR, CL, and RL
  ZRL = ZR + ZL ;
  
  CT = (1/w^2-(LT+imag(w*M^2/ZRL))*real(CPT))/(LT+imag(w*M^2/ZRL)-real(CPT)*((RT+real(w^2*M^2/ZRL))^2+(w*LT+imag(w^2*M^2/ZRL))^2))-real(CPT);
  if CT > 0
    sCT = s*CT;
    ZCT = 1/sCT;
  else
    ZCT = 0;
  end
  
  sCPT = s*CPT;         %Admittances of capacitances at frequency f

        
  ZT = RT + sLT - sM;                                       %ZT is calculated by adding impedances of RT and LT
                    
  VM = VR / ZL * ( ZR + ZL );                               %VM is calculated from VR using voltage divider
  VTin = VM + ZT * (( VM / sM ) + ( VM / ( ZR + ZL ) ) );   %VTin is calculated from VM + voltage drop across ZT
  VT = VTin + ZCT * (VTin*sCPT + (VTin-VM)/ZT);           %VT is calculated by VTin + voltage drop across CT
        
  IT = VTin * sCPT + ( VTin - VM ) / ZT;                    %IT is calculated by adding current through sM and current through (ZR+ZL)
        
  PT = sqrt( VT * IT * conj( VT * IT ) );                   %Input power is defined as the magnitude of the complex power on the input
  PL = sqrt( VR * VR / RL * conj( VR * VR / RL ) );         %Output power is defined as the magnitude of the complex power on the load resistance
  
  eta = PL / PT;                                            %Power transfer efficiency is the ratio of the input and output power
  
end

function pte_plot3(x,pte,xlabel_,ylabel_,xmin,xmax,x_res,ymin,ymax,y_res)
  grid on;
  box off
  fs=9;lw=0.5;
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  plot(x,pte,'ko-','linewidth',lw,'markersize',3,'markerfacecolor','none')
 
  if(abs(1-exist('xmin'))>0)
    xmin=0;
  end
  if(abs(1-exist('xmax'))>0 || xmax=='a')
    xmax=max(xlim);
    xlim([xmin xmax])
  else
    xlim([xmin xmax])
  end

  if(abs(1-exist('ymin'))>0)
    if min(k) > 0
      ymin=0;
    else
      ymin=min(k);
    end
  elseif(ymin=='a')
    ymin=min(k);
  end
  if(abs(1-exist('ymax'))>0 || ymax=='a')
    ymax=max(ylim);
    ylim([ymin ymax])  
  else
    ylim([ymin ymax])  
  end

  if(abs(1-exist('x_res'))>0) || x_res=='a'
    x_res=5;
  else
    set(gca,'xtick',[xmin:(xmax-xmin)/x_res:xmax]);    
  end 
  if(abs(1-exist('y_res'))>0 || y_res=='a')
    y_res=3;
  else
    set(gca,'ytick',[ymin:(ymax-ymin)/y_res:ymax]);    
  end
  
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');
  aa=get(gca,'position');
  set(gca,'ticklength',[0.001*1/aa(3) 0]); 
  set(gca,'tickdir','out');
  x = get(gca,'position');
  x(3)=x(3);
  set(gca,'position',x);
  
  set(gca,'xminorgrid','off','yminorgrid','off','minorgridlinestyle','-');
  major_gc_=0;
  minor_gc_=0.75;
  majorgridcolor_ = [major_gc_ major_gc_ major_gc_];
  minorgridcolor_ = [minor_gc_ minor_gc_ minor_gc_];
  set(gca,'gridcolor', majorgridcolor_, 'minorgridcolor',minorgridcolor_);
  xlabel(xlabel_);ylabel(ylabel_);
  set(gca,'fontsize',fs,'linewidth',lw,'fontname','arial');  
  grid on;
end

