clear;
close all;
clc;

% Fixed values
f=13.56e6;
RL=100;
VR=5;
z=15000;

% Coil parameters
odT = 52000;
idT = 10000; 
sT = 1000;
tT = 17.5;
nT = 9;
sigmaT = 58000000;
RsheetT = 1/sigmaT/(tT*1e-6);

odR = 5000;
idR = 3750;
sR = 5;
tR = 2.05;
nR = 7;
sigmaR = 38000000;
RsheetR = 1/sigmaR/(tR*1e-6);

sR = linspace(1,60,60)
for i=1:60
    RT(i) = Res_prox(sT, nT, odT, idT, RsheetT, RsheetT, 0, tT, f);
    [LT(i), wT(i)] = wheelerTurns(sT, nT, odT, idT);
    RR(i) = Res(sR(i), nR, odR, idR, RsheetR, RsheetR, 0, tT, f);
    [LR(i), wR(i)] = wheelerTurns(sR(i), nR, odR, idR);
    k(i) = kFilament(sR(i), nR, odR, idR, sT, nT, odT, idT, z);
    [PTE_matlab(i),CT_matlab(i),CL_matlab(i),VT_matlab(i),VM_matlab(i),IT_matlab(i),IR_matlab(i)] = maxPTE_Csweep(RT(i),LT(i),1e-12,10e-12,100e-12,RR(i),LR(i),1e-12,10e-12,200e-12,k(i),f,RL,VR,5);
end


plot(sR,PTE_matlab);

xlabel('rx spacing')
ylabel('\eta')
ylim([0 max(ylim)]);
cd "P:/Bioelectronics-Kim/allusers/Paper/22-2020-IEEE_TBioCAS_antennadesign/Figure/";
cd "Figure 4 - k and PTE vs design parameters/Matlab_Results/";
saveas(gca,"eta_vs_sr_matlab.png")

figure('position',[100 100 800 600])
fs=20;

subplot(3,4,1)
txn = dlmread('P:\BioElectronics-Kim\allusers\Paper\22-2020-IEEE_TBioCAS_antennadesign\Figure\Figure 2 - k decoupled from other parameters\HFSS_Results\txn.csv',",");      
sim_n = [5, 26];
kplot_( txn(sim_n(1):sim_n(2),14), txn(sim_n(1):sim_n(2),18));
text(-0,0,"0 ",'horizontalalignment','right','verticalalignment','bottom','fontsize',fs)
text(-0,2,"2 ",'horizontalalignment','right','verticalalignment','top','fontsize',fs)
text(-6,1,"Tx coil{\n}k, %",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'rotation',90)



cd "P:/BioElectronics-Kim/allusers/Paper/22-2020-IEEE_TBioCAS_antennadesign/Figure/";
cd "Figure 4 - k and PTE vs design parameters/Matlab_Results/";
saveas(gca,strcat("eta_subplots.png"))