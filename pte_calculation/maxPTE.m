% Input: effective R and L of the TX and RX coils, tuning cap sweep range and resolution, k between TX and RX coils, frequency, load resistance, and target received voltage.
% Output: maximum power transfer efficiency, tuning caps, voltages and currents of TX and RX ports
%
%
% Circuit model
%
%     VT               VM                   VR
%     *-||---^^^---mmm-*-mmm---^^^-----------*
%       CT   RT   LT-M | LR-M   RR  |        |
%                      3           ---      <
%                    M 3        CL ---   RL <   
%                      |            |        |
%     ----------------------------------------
%
% 
% The circuit will be simplified as below before calculating power transfer efficiency
%
%     VT               VM                   VR
%     *-----|ZT|-------*-------|ZR|---------*
%     -->              |     -->             |       
%     IT               |     IR              |
%                     |sM|                  |ZL|
%                      |                     | 
%     ----------------------------------------
%

function [eta_max,CT_,CL_,VT_,VM_,IT_,IR_] = maxPTE(RT,LT,CT_low,CT_res,CT_high,RR,LR,CL_low,CL_res,CL_high,k,f,RL,VR)
M = k * sqrt( LT * LR );
s = 1i * 2 * pi * f; sLT = s * LT; sLR = s * LR; sM = s * M;      %Impedances of inductances at frequency f    
          
eta_max=0;
for CT = CT_low : CT_res : CT_high                                %CT sweep
    for CL = CL_low : CL_res : CL_high                            %CL sweep
        sCT = s * CT; sCL = s * CL;                               %Admittances of capacitances at frequency f
        
        ZT = RT + sLT - sM + 1 / sCT;                             %ZT is calculated by adding impedances of RT, LT, and CT
        ZR = RR + sLR - sM;                                       %ZR is calculated by adding impedances of RR and LR
        ZL = 1 / ( sCL + 1 / RL );                                %ZL is calculated by adding admittances of CL and RL
                    
        VM = VR / ZL * ( ZR + ZL );                               %VM is calculated from VR using voltage divider
        VT = VM * ( 1+ ( ZT / sM ) + ( ZT / ( ZR + ZL ) ) );      %VT is calculated from VM using mesh equation
 
        IT = VM / sM + VM / (ZR + ZL);                            %IT is calculated by adding current through sM and current through (ZR+ZL)
        IR = VM / ( ZR + ZL );                                    %IR is calculated as dividing VM by the equvalent impedance of ZR and ZL
        
        PT = sqrt( VT * IT * conj( VT * IT ) );                   %Input power is defined as the magnitude of the complex power on the input
        PL = sqrt( VR * VR / RL * conj( VR * VR / RL ) );         %Output power is defined as the magnitude of the complex power on the load resistance
 
        eta = PL / PT;                                            %Power transfer efficiency is the ratio of the input and output power
                                       
        if eta > eta_max
             eta_max = eta; CT_ = CT; CL_ = CL; VT_ = VT; VM_ = VM; IT_ = IT; IR_ = IR; %Find the maximum power transfer efficiency and save the outputs
        end
    end
end
endfunction