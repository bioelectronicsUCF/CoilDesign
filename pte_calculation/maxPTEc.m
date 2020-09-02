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

function [eta,CT,CL,VT,VM,IT] = maxPTEc_sweep_cap(RT,LT,CPT,RR,LR,CPR,M,f,RL,VR)

s = 1i * 2 * pi * f; sLT = s * LT; sLR = s * LR; sM = s * M;      %Impedances of inductances at frequency f    
w = 2 * pi * f;
          
CT = real(-1/s / (1/(s*CPT+1/(RT+sLT))));
CL_max = 1/(RR^2/LR+w^2*LR)
CL_init = real(CPR)
if CL_max > CL_init
  CL = CL_max - CL_init;
else
  CL = 0;
end

  sCT = s*CT; sCPT = s*CPT; sCR = s*(CL+CPR);        %Admittances of capacitances at frequency f

        
  ZT = RT + sLT - sM;                                       %ZT is calculated by adding impedances of RT and LT
  ZR = RR + sLR - sM;                                       %ZR is calculated by adding impedances of RR and LR
  ZL = 1 / ( sCR + 1 / RL );                         %ZL is calculated by adding admittances of CPR, CL, and RL
                    
  VM = VR / ZL * ( ZR + ZL );                               %VM is calculated from VR using voltage divider
  VTin = VM + ZT * (( VM / sM ) + ( VM / ( ZR + ZL ) ) );   %VTin is calculated from VM + voltage drop across ZT
  VT = VTin + 1/sCT * (VTin*sCPT + (VTin-VM)/ZT);           %VT is calculated by VTin + voltage drop across CT
        
  IT = VTin * sCPT + ( VTin - VM ) / ZT;                    %IT is calculated by adding current through sM and current through (ZR+ZL)
        
  PT = sqrt( VT * IT * conj( VT * IT ) );                   %Input power is defined as the magnitude of the complex power on the input
  PL = sqrt( VR * VR / RL * conj( VR * VR / RL ) );         %Output power is defined as the magnitude of the complex power on the load resistance
 
  eta = PL / PT;                                            %Power transfer efficiency is the ratio of the input and output power
  
endfunction