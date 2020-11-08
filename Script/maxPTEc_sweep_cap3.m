
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
