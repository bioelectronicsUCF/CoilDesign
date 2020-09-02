function PTE = ptecal_sr_sweep(p, s)
% p: parameters, s: sweep info
% Fixed values, p. f, RL, VR, z
% Coil parameters, p. odT, idT, sT, tT, nT, sigmaT, odR,idR,sR,tR,nR,sigmaR
p.RsheetT = 1/p.sigmaT/(p.tT*1e-6);
p.RsheetR = 1/p.sigmaR/(p.tR*1e-6);

p.sR = linspace(s.start_,s.end_,s.resolution_)
it_res = 5 %number of iterations
for i=1:s.resolution_
    RT(i) = Res_prox(p.sT, p.nT, p.odT, p.idT, p.RsheetT, p.RsheetT, 0, p.tT, p.f);
    [LT(i), wT(i)] = wheelerTurns(p.sT, p.nT, p.odT, p.idT);
    RR(i) = Res(p.sR(i), p.nR, p.odR, p.idR, p.RsheetR, p.RsheetR, 0, p.tT, p.f);
    [LR(i), wR(i)] = wheelerTurns(p.sR(i), p.nR, p.odR, p.idR);
    k(i) = kFilament(p.sR(i), p.nR, p.odR, p.idR, p.sT, p.nT, p.odT, p.idT, p.z);
    [PTE(i),CT(i),CL(i),VT(i),VM(i),IT(i),IR(i)] = maxPTE_Csweep(RT(i),LT(i),1e-12,10e-12,100e-12,RR(i),LR(i),1e-12,10e-12,200e-12,k(i),p.f,p.RL,p.VR,it_res);
end


plot(sR,PTE_matlab);

xlabel('rx spacing')
ylabel('\eta')
ylim([0 max(ylim)]);
cd "P:/Bioelectronics-Kim/allusers/Paper/22-2020-IEEE_TBioCAS_antennadesign/Figure/";
cd "Figure 4 - k and PTE vs design parameters/Matlab_Results/";
saveas(gca,"eta_vs_sr_matlab.png")
