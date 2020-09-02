function C = Initialize7()

% Fixed values
  C.f=13.56e6;
  C.RL=100;
  C.VR=5;
  C.z=15000;

% Coil parameters
  C.tod = 50000;
  C.tid = 10000; 
  C.ts = 500;
  C.tt = 500;
  C.tn = 20;
  C.tsigma = 58000000;
  C.tRsheet = 1/C.tsigma/(C.tt*1e-6);
  C.tepsilon_t_sub = 0;
  C.tepsilon_t_mm = 3.5*(1-0.008j);
  C.tepsilon_t_mmf = 3.5*(1-0.008j);
  C.td = 1000; 
  C.th = 1000;
  C.tcon = 20000;

  C.rod = 15000;
  C.rid = 13000;
  C.rs = 100;
  C.rt = 30;
  C.rh = 10;
  C.rn = 3;
  C.rsigma = 37700000;
  C.rRsheet = 1/C.rsigma/(C.rt*1e-6);
  C.repsilon_r_sub = 4.1;
  C.rt_ins = 10;
  C.repsilon_r_mm = 4.1;
  C.repsilon_r_mmf = 4.1;
  C.rd = 5;

end