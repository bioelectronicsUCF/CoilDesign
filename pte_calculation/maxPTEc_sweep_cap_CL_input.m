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

function [eta,CT,CL,VT,VM,IT] = maxPTEc_sweep_cap_CL_input(RT,LT,CPT,RR,LR,CPR,M,f,RL,VR,CL)

s = 1i * 2 * pi * f; sLT = s * LT; sLR = s * LR; sM = s * M;      %Impedances of inductances at frequency f    
w = 2 * pi * f;
          
%CT = real(-1/s / (1/(s*CPT+1/(RT+sLT))));

  
  sCR = s*(CL+CPR);
  
  ZR = RR + sLR - sM;                                       %ZR is calculated by adding impedances of RR and LR
  ZL = 1 / ( sCR + 1 / RL );                                %ZL is calculated by adding admittances of CPR, CL, and RL
  ZRL = ZR + ZL + sM;
  
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