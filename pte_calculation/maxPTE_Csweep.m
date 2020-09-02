% incrementally enhance resolution of the sweep 

function [eta_max,CT,CL,VT,VM,IT,IR] = maxPTE_Csweep(RT,LT,CT_low,CT_res,CT_high,RR,LR,CL_low,CL_res,CL_high,k,f,RL,VR,sweep_n)

for i = 1 : sweep_n
  [eta_max,CT,CL,VT,VM,IT,IR] = maxPTE(RT,LT,CT_low,CT_res,CT_high,RR,LR,CL_low,CL_res,CL_high,k,f,RL,VR);
  CT_low = CT - CT_res;
  CT_high = CT + CT_res;
  CT_res = CT_res / 10;
  CL_low = CL - CL_res;
  CL_high = CL + CL_res;
  CL_res = CL_res / 10;
end

endfunction