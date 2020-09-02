function eta_max=PTE_from_Z_Csweep(z11,z12,z22)
z11;z12;z22;
CT_start = 1e-15;
CT_end = 200000e-12;
CT_inc = 1000e-12;

CL_start = 1e-12;
CL_end = 200000e-12;
CL_inc = 1000e-12;

eta_max=0;CT_max=0;CL_max=0;VT_max=0;VM_max=0;IT_max=0;IR_max=0;
for iteration = 1 : 9
    [eta_max,CT_max,CL_max,VT_max,VM_max,IT_max,IR_max] = PTE_from_Z(z11,z12,z22,CL_start,CL_end,CL_inc,CT_start,CT_end,CT_inc);
    CT_start = max(CT_start,CT_max-CT_inc);
    CT_end   = min(CT_end,CT_max+CT_inc) ;
    CT_inc   = CT_inc/10;
    CL_start = max(CL_start,CL_max-CL_inc);
    CL_end   = min(CL_end,CL_max+CL_inc) ;
    CL_inc   = CL_inc/10;
end
ans='aa'
end