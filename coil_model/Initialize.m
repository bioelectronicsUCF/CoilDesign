function [f, RL, VR, z, C] = Initialize();

% Fixed values
  f=13.56e6;
  RL=100;
  VR=5;
  z=15000;

% Coil parameters
  C.tod = 50000;
  C.tid = 10000; 
  C.ts = 1000;
  C.tt = 175;
  C.tn = 10;
  C.tsigma = 58000000;
  C.tRsheet = 1/C.tsigma/(C.tt*1e-6);
  C.tepsilon_t_sub = 0;
  C.tepsilon_t_mm = 3.5*(1-0.008j);
  C.tepsilon_t_mmf = 3.5*(1-0.008j);
  C.td = 1000; 
  C.th = 1000;
  C.tcon = 20000;

  C.rod = 5000;
  C.rid = 4000;
  C.rs = 5;
  C.rt = 0.925;
  C.rh = 6.155;
  C.rn = 8;
  C.rRsheet = 0.051;
  C.rsigma = 1/C.rRsheet /C.rt*1e6;
  C.rt_insulator1 = 1.39;    % in um (field oxide thickness)
  C.repsilon_r_insulator1 = 3.9; % relative permittivity of field oxide
  C.rt_insulator2 = 4.765;   % in um (ILD thickness)
  C.repsilon_r_insulator2 = 4.1; % relative permittivity of ILD
  C.rt_insulator3 = 0.54;
  C.repsilon_r_insulator3 = 7.9;
  C.repsilon_r_sub = (1.39*3.9+4.765*4.1+0.54*7.9)/(1.39+4.765+0.54);
  C.rt_ins = 1.39+4.765+0.54;
  C.repsilon_r_mm = C.repsilon_r_sub;
  C.repsilon_r_mmf = C.repsilon_r_sub;
  %C.repsilon_r_sub = C.repsilon_r_insulator1 * C.rt_insulator1/C.rh + C.repsilon_r_insulator2 * C.rt_insulator2/C.rh;
  %C.repsilon_r_mm = C.repsilon_r_insulator2;
  %C.repsilon_r_mmf = 1*0.5;
  C.rd = 4.28;

endfunction