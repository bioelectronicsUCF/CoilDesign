function [RT_turns, LT_turns, Cp_turns] = Z_tot(C,f,csub,gsub)

RT_turns = Res_tot_turns(C.ts, C.tn, C.tod, C.tid, C.tRsheet, C.tRsheet, C.tcon, C.tt, f);
LT_turns = Greenhouse_turns(C.ts, C.tn, C.tod, C.tid, C.tt);
Cp_turns = Capp_turns(C.ts, C.tn, C.tod, C.tid, C.tepsilon_t_mmf, C.tt);

endfunction