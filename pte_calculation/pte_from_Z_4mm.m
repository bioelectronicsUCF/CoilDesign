clear
clc;
close all;

%import_=importdata('P:\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\Data\2ozsubfrz4idtxs_wire.csv');
import_=importdata('..\..\HFSS_Results\4id_dense.csv');
%import_=importdata('P:\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\Data\2ozsubfrz4idtxs.csv');
%import_=importdata('P:\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\Data\2ozsubfrz4id.csv');
%import_=importdata('P:\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\Data\2ozsubfrz.csv');
%import_=importdata('P:\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\Data\1ozsubfrz.csv');
import_size=size(import_.data);
m_sweep=import_size(1);
n_sweep=(import_size(2)-1)/8;

z11=import_.data(1:m_sweep,2:n_sweep+1)+1i*import_.data(1:m_sweep,n_sweep+2:2*n_sweep+1);
z12=import_.data(1:m_sweep,2*n_sweep+2:3*n_sweep+1)+1i*import_.data(1:m_sweep,3*n_sweep+2:4*n_sweep+1);
z21=import_.data(1:m_sweep,4*n_sweep+2:5*n_sweep+1)+1i*import_.data(1:m_sweep,5*n_sweep+2:6*n_sweep+1);
z22=2*import_.data(1:m_sweep,6*n_sweep+2:7*n_sweep+1)+1i*import_.data(1:m_sweep,7*n_sweep+2:8*n_sweep+1);

f=13.56e6;
w=2*pi*f;
s=1i*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CT      Z11-Z12  VM   Z22-Z12
%    -------||------[:::::]-------[:::::]------------------
%   +   -->                   |     -->       |    -->    | +
%       Iin                   -     I2        |    IL     |
%                            | |             ---          <
%  Vin                   Z12 | |          CL ---       RL < Vout
%                             -               |           |
%   -                         |               |           | -
%    ------------------------------------------------------


Vout = 5;
RL = 100;
IL = Vout/RL;

Pin_min=ones(m_sweep, n_sweep)*1e6;
CL_max=zeros(m_sweep, n_sweep);
CT_max=zeros(m_sweep, n_sweep);

for CL=50e-12:1e-12:300e-12
    
    I2=IL+Vout*s*CL;
    VM=Vout+I2*(z22-z12);
    Iin=I2+VM./z12;
    CL;
    for CT=20e-12:0.1e-12:60e-12
    
        Vin=VM+Iin.*(1/(s*CT)+z11-z12);
        Pin=abs(Iin.*Vin);
        
        for m=1:1:m_sweep
            for n=1:1:n_sweep
                if Pin(m,n)<Pin_min(m,n)
                    Pin_min(m,n)=Pin(m,n);
                    CL_max(m,n)=CL;
                    CT_max(m,n)=CT;
                end
            end
        end
    end
end

I2=IL+Vout*s.*CL_max;
VM=Vout+I2.*(z22-z12);
Iin=I2+VM./z12;
Vin=VM+Iin.*(1./(s.*CT_max)+z11-z12);
Pin=abs(Iin.*Vin);

surf(9:13,6:5+m_sweep,0.25./Pin)
xlabel('Tx Turns');
ylabel('Rx Turns');
title('PTE to 100 ohm load vs Tx and Rx turns');
grid on;

figure()
hold on;

eta = 0.25./Pin
txn=linspace(9,8+n_sweep,n_sweep);
rxn=linspace(6,5+m_sweep,m_sweep);

for txn_index = 1 : 5
  
  [max_pte max_pte_indice] = max(eta(:,txn_index));
  [min_pte min_pte_indice] = min(eta(:,txn_index));
  
  red_point = min_pte + (max_pte-min_pte)/3;
  green_point = min_pte + (max_pte-min_pte)/3*2;
 
  for rxn_index = 1:m_sweep
    if (eta(txn_index,rxn_index)>green_point)
      plot(txn(txn_index),rxn(rxn_index),'go',"markersize",15)
    elseif (eta(txn_index,rxn_index)>red_point)
      plot(txn(txn_index),rxn(rxn_index),'yo',"markersize",15)
    else
      plot(txn(txn_index),rxn(rxn_index),'ro',"markersize",15)
    endif
  endfor
  plot(txn(txn_index),rxn(max_pte_indice),'kx',"markersize",15);
  x(txn_index)=txn(txn_index);
  y(txn_index)=rxn(max_pte_indice);
endfor
