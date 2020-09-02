function [eta_max,CT_max,CL_max,VT_max,VM_max,IT_max,IR_max]=PTE_from_Z(z11,z12,z22,CL_start,CL_end,CL_inc,CT_start,CT_end,CT_inc)

f=13.56e6;
w=2*pi*f;
s=1i*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CT      Z11-Z12  VM   Z22-Z12
%    -------||------[:::::]-------[:::::]------------------
%   +   -->                   |     -->       |    -->    | +
%       Iin                   -     I2        |    Iout   |
%                            | |             ---          <
%  Vin                   Z12 | |          CL ---       RL < Vout
%                             -               |           |
%   -                         |               |           | -
%    ------------------------------------------------------


Vout = 5;
RL   = 100;
Iout = Vout/RL;

eta_max = 0;
CL_max  = 0;
CT_max  = 0;

#CT = abs(1/(w*imag(z11)));
for CT = CT_start:CT_inc:CT_end
for CL = CL_start:CL_inc:CL_end
    
  I2  = Iout+Vout*s*CL;
  VM  = Vout+I2*(z22-z12);
  Iin = I2+VM/z12;
    
  Vin = VM+Iin*(1/(s*CT)+z11-z12);
  Pin = abs(Iin*Vin);
  eta = (abs(Iout*Vout))/Pin;
        
  if(eta>eta_max)
            
    eta_max = eta;
    CT_max = CT;
    CL_max = CL;
    VT_max = Vin;
    VM_max = VM;
    IT_max = Iin;
    IR_max = Iout;
            
  end
end
end
endfunction