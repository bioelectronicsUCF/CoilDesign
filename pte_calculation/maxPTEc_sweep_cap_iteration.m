% incrementally enhance resolution of the sweep 

function [eta_max,CT,CL,VT,VM,IT,IR] = maxPTEc_sweep_cap_iteration(RT,LT,CPT,CT,RR,LR,CPR,CL,M,f,RL,VR)

CT_start = 1e-15;
CT_end = 2000e-12;
CT_inc = 10e-12;

CL_start = 1e-12;
CL_end = 2000e-12;
CL_inc = 10e-12;

for i = 1 : 9
  [eta_max,CT,CL,VT,VM,IT,IR] = maxPTEc_sweep_cap2(RT,LT,CPT,CT,RR,LR,CPR,CL,M,f,RL,VR,CL_start,CL_end,CL_inc,CT_start,CT_end,CT_inc);
  CT_low = CT - CT_res;
  CT_high = CT + CT_res;
  CT_res = CT_res / 10;
  CL_low = CL - CL_res;
  CL_high = CL + CL_res;
  CL_res = CL_res / 10;
end

endfunction