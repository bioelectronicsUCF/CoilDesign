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

