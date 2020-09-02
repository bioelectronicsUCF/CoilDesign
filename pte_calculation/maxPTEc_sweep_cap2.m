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

function [eta_max,CT_max,CL_max,VT_max,VM_max,IT_max] = maxPTEc_sweep_cap2(RT,LT,CPT,CT,RR,LR,CPR,CL,M,f,RL,VR,CL_start,CL_end,CL_inc,CT_start,CT_end,CT_inc)

s = 1i * 2 * pi * f; sLT = s * LT; sLR = s * LR; sM = s * M;      %Impedances of inductances at frequency f    
w = 2 * pi * f;

eta_max = 0;
CL_max  = 0;
CT_max  = 0;

for CT = CT_start:CT_inc:CT_end
for CL = CL_start:CL_inc:CL_end          
  sCT = s*CT; sCPT = s*CPT; sCR = s*(CL+CPR);               %Admittances of capacitances at frequency f
        
  ZT = RT + sLT - sM;                                       %ZT is calculated by adding impedances of RT and LT
  ZR = RR + sLR - sM;                                       %ZR is calculated by adding impedances of RR and LR
  ZL = 1 / ( sCR + 1 / RL );                                %ZL is calculated by adding admittances of CPR, CL, and RL
                    
  VM = VR / ZL * ( ZR + ZL );                               %VM is calculated from VR using voltage divider
  VTin = VM + ZT * (( VM / sM ) + ( VM / ( ZR + ZL ) ) );   %VTin is calculated from VM + voltage drop across ZT
  VT = VTin + 1/sCT * (VTin*sCPT + (VTin-VM)/ZT);           %VT is calculated by VTin + voltage drop across CT
        
  IT = VTin * sCPT + ( VTin - VM ) / ZT;                    %IT is calculated by adding current through sM and current through (ZR+ZL)
        
  PT = sqrt( VT * IT * conj( VT * IT ) );                   %Input power is defined as the magnitude of the complex power on the input
  PL = sqrt( VR * VR / RL * conj( VR * VR / RL ) );         %Output power is defined as the magnitude of the complex power on the load resistance
 
  eta = PL / PT;                                            %Power transfer efficiency is the ratio of the input and output power
        
  if(eta>eta_max)
    eta_max = eta;
    CT_max = CT;
    CL_max = CL;
    VT_max = VT;
    VM_max = VM;
    IT_max = Iin;
    IR_max = Iout;
  end  
end
end
endfunction