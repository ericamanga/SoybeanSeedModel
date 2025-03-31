% Setting different biomass composition according to the conditions each
% plant was submitted in the field:

% Define paths for model reconstruction.
clear; % clean workspace
clc; % clean command window
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
 
%% Load curated model 'Glycine max' GmaxGEM with RAVEN toolbox
% The curation was obtained with manualCuration.m script

% Load the curated model
load([data '/mat/soybeanSeedModel.mat'],'model')

% Firsly, we set all biomass composition to either 0 or fixed as
% experimental data provided. First setting all biomass to zero:
BiomassRxns = model.rxns(endsWith(model.rxns,'_biomass'));
model = setParam(model, 'eq', BiomassRxns,0);

% Then, setting the different biomass composition:

%% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
model = setParam(model,'eq','Starch_biomass',-0.000410487);
model = setParam(model,'eq','sGLC_biomass',-2.57828E-05);
model = setParam(model,'eq','sFRU_biomass',-4.12566E-05);
model = setParam(model,'eq','sSUCROSE_biomass',-0.000142125);
model = setParam(model,'eq','Raffinose_biomass',-3.94316E-06);
model = setParam(model,'eq','Stachyose_biomass',-5.53559E-06);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-6.03772E-08);
model = setParam(model,'eq','Palmitate_biomass',-3.96094E-05);
model = setParam(model,'eq','Stearate_biomass',-9.41741E-06);
model = setParam(model,'eq','Oleate_biomass',-4.11526E-05);
model = setParam(model,'eq','Linoleate_biomass',-9.50984E-05);
model = setParam(model,'eq','Linolenate_biomass',-3.90439E-05);
model = setParam(model,'eq','Arachidate_biomass',-7.79281E-07);
model = setParam(model,'eq','Glycerol_biomass',-7.50538E-05);
model = setParam(model,'eq','CELLWALL_biomass',-0.001421555);
model = setParam(model,'eq','pHIS_biomass',-2.77399E-05);
model = setParam(model,'eq','pILE_biomass',-6.66464E-05);
model = setParam(model,'eq','pLEU_biomass',-0.000130934);
model = setParam(model,'eq','pLYS_biomass',-9.37275E-05);
model = setParam(model,'eq','pMET_biomass',-1.93309E-05);
model = setParam(model,'eq','pCYS_biomass',-4.33523E-05);
model = setParam(model,'eq','pPHE_biomass',-8.29882E-05);
model = setParam(model,'eq','pTYR_biomass',-3.80394E-05);
model = setParam(model,'eq','pTHR_biomass',-7.39308E-05);
model = setParam(model,'eq','pTRP_biomass',-1.7563E-05);
model = setParam(model,'eq','pVAL_biomass',-8.28074E-05);
model = setParam(model,'eq','pARG_biomass',-9.23019E-05);
model = setParam(model,'eq','pALA_biomass',-0.000112666);
model = setParam(model,'eq','pASP_biomass',-0.00011482);
model = setParam(model,'eq','pASN_biomass',-0.000115812);
model = setParam(model,'eq','pGLU_biomass',-0.000135875);
model = setParam(model,'eq','pGLN_biomass',-0.000156005);
model = setParam(model,'eq','pGLY_biomass',-0.000125153);
model = setParam(model,'eq','pPRO_biomass',-0.000160806);
model = setParam(model,'eq','pSER_biomass',-0.000115727);

save([data '/mat/wr1_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.001153073);
model = setParam(model,'eq','sGLC_biomass',-6.17857E-05);
model = setParam(model,'eq','sFRU_biomass',-9.66636E-05);
model = setParam(model,'eq','sSUCROSE_biomass',-0.00044332);
model = setParam(model,'eq','Raffinose_biomass',-2.54231E-06);
model = setParam(model,'eq','Stachyose_biomass',0);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-1.71105E-07);
model = setParam(model,'eq','Palmitate_biomass',-0.000149469);
model = setParam(model,'eq','Stearate_biomass',-4.26232E-05);
model = setParam(model,'eq','Oleate_biomass',-0.00020708);
model = setParam(model,'eq','Linoleate_biomass',-0.000358365);
model = setParam(model,'eq','Linolenate_biomass',-0.000102964);
model = setParam(model,'eq','Arachidate_biomass',-3.28809E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000287987);
model = setParam(model,'eq','CELLWALL_biomass',-0.003740036);
model = setParam(model,'eq','pHIS_biomass',-7.68782E-05);
model = setParam(model,'eq','pILE_biomass',-0.000184703);
model = setParam(model,'eq','pLEU_biomass',-0.000362869);
model = setParam(model,'eq','pLYS_biomass',-0.000259756);
model = setParam(model,'eq','pMET_biomass',-5.35737E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000120146);
model = setParam(model,'eq','pPHE_biomass',-0.000229993);
model = setParam(model,'eq','pTYR_biomass',-0.000105422);
model = setParam(model,'eq','pTHR_biomass',-0.000204891);
model = setParam(model,'eq','pTRP_biomass',-4.86739E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000229492);
model = setParam(model,'eq','pARG_biomass',-0.000255805);
model = setParam(model,'eq','pALA_biomass',-0.000312243);
model = setParam(model,'eq','pASP_biomass',-0.000318213);
model = setParam(model,'eq','pASN_biomass',-0.000320959);
model = setParam(model,'eq','pGLU_biomass',-0.000376563);
model = setParam(model,'eq','pGLN_biomass',-0.00043235);
model = setParam(model,'eq','pGLY_biomass',-0.000346849);
model = setParam(model,'eq','pPRO_biomass',-0.000445655);
model = setParam(model,'eq','pSER_biomass',-0.000320726);


save([data '/mat/wr2_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.002208058);
model = setParam(model,'eq','sGLC_biomass',-0.000114932);
model = setParam(model,'eq','sFRU_biomass',-0.000151478);
model = setParam(model,'eq','sSUCROSE_biomass',-0.000735542);
model = setParam(model,'eq','Raffinose_biomass',-1.13587E-05);
model = setParam(model,'eq','Stachyose_biomass',0);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-3.07137E-07);
model = setParam(model,'eq','Palmitate_biomass',-0.000221672);
model = setParam(model,'eq','Stearate_biomass',-6.80679E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000364437);
model = setParam(model,'eq','Linoleate_biomass',-0.00043201);
model = setParam(model,'eq','Linolenate_biomass',-8.85133E-05);
model = setParam(model,'eq','Arachidate_biomass',-5.71899E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000393576);
model = setParam(model,'eq','CELLWALL_biomass',-0.007448494);
model = setParam(model,'eq','pHIS_biomass',-0.000136724);
model = setParam(model,'eq','pILE_biomass',-0.000328484);
model = setParam(model,'eq','pLEU_biomass',-0.000645341);
model = setParam(model,'eq','pLYS_biomass',-0.000461961);
model = setParam(model,'eq','pMET_biomass',-9.52777E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000213673);
model = setParam(model,'eq','pPHE_biomass',-0.000409029);
model = setParam(model,'eq','pTYR_biomass',-0.000187487);
model = setParam(model,'eq','pTHR_biomass',-0.000364387);
model = setParam(model,'eq','pTRP_biomass',-8.65638E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000408138);
model = setParam(model,'eq','pARG_biomass',-0.000454934);
model = setParam(model,'eq','pALA_biomass',-0.000555306);
model = setParam(model,'eq','pASP_biomass',-0.000565923);
model = setParam(model,'eq','pASN_biomass',-0.000570808);
model = setParam(model,'eq','pGLU_biomass',-0.000669696);
model = setParam(model,'eq','pGLN_biomass',-0.00076891);
model = setParam(model,'eq','pGLY_biomass',-0.000616852);
model = setParam(model,'eq','pPRO_biomass',-0.000792573);
model = setParam(model,'eq','pSER_biomass',-0.000570393);

save([data '/mat/wr3_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.001881329);
model = setParam(model,'eq','sGLC_biomass',-0.000128833);
model = setParam(model,'eq','sFRU_biomass',-0.000122015);
model = setParam(model,'eq','sSUCROSE_biomass',-0.000690668);
model = setParam(model,'eq','Raffinose_biomass',-1.02338E-05);
model = setParam(model,'eq','Stachyose_biomass',0);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-3.22931E-07);
model = setParam(model,'eq','Palmitate_biomass',-0.000252694);
model = setParam(model,'eq','Stearate_biomass',-8.26977E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000459744);
model = setParam(model,'eq','Linoleate_biomass',-0.000536851);
model = setParam(model,'eq','Linolenate_biomass',-9.67426E-05);
model = setParam(model,'eq','Arachidate_biomass',-7.28451E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000478779);
model = setParam(model,'eq','CELLWALL_biomass',-0.007889106);
model = setParam(model,'eq','pHIS_biomass',-0.00014925);
model = setParam(model,'eq','pILE_biomass',-0.000358581);
model = setParam(model,'eq','pLEU_biomass',-0.000704469);
model = setParam(model,'eq','pLYS_biomass',-0.000504287);
model = setParam(model,'eq','pMET_biomass',-0.000104007);
model = setParam(model,'eq','pCYS_biomass',-0.00023325);
model = setParam(model,'eq','pPHE_biomass',-0.000446505);
model = setParam(model,'eq','pTYR_biomass',-0.000204666);
model = setParam(model,'eq','pTHR_biomass',-0.000397773);
model = setParam(model,'eq','pTRP_biomass',-9.4495E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000445533);
model = setParam(model,'eq','pARG_biomass',-0.000496617);
model = setParam(model,'eq','pALA_biomass',-0.000606184);
model = setParam(model,'eq','pASP_biomass',-0.000617774);
model = setParam(model,'eq','pASN_biomass',-0.000623107);
model = setParam(model,'eq','pGLU_biomass',-0.000731055);
model = setParam(model,'eq','pGLN_biomass',-0.00083936);
model = setParam(model,'eq','pGLY_biomass',-0.000673369);
model = setParam(model,'eq','pPRO_biomass',-0.00086519);
model = setParam(model,'eq','pSER_biomass',-0.000622654);

save([data '/mat/wr4_c1-c2.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.001160589);
model = setParam(model,'eq','sGLC_biomass',-3.67885E-05);
model = setParam(model,'eq','sFRU_biomass',-1.07489E-05);
model = setParam(model,'eq','sSUCROSE_biomass',-0.000474884);
model = setParam(model,'eq','Raffinose_biomass',3.30734E-06);
model = setParam(model,'eq','Stachyose_biomass',4.42847E-06);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-1.9923E-07);
model = setParam(model,'eq','Palmitate_biomass',-0.000168646);
model = setParam(model,'eq','Stearate_biomass',-7.25057E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000330521);
model = setParam(model,'eq','Linoleate_biomass',-0.000510342);
model = setParam(model,'eq','Linolenate_biomass',-0.000104181);
model = setParam(model,'eq','Arachidate_biomass',-5.11927E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000397171);
model = setParam(model,'eq','CELLWALL_biomass',-0.004386975);
model = setParam(model,'eq','pHIS_biomass',-9.14951E-05);
model = setParam(model,'eq','pILE_biomass',-0.000219821);
model = setParam(model,'eq','pLEU_biomass',-0.000431861);
model = setParam(model,'eq','pLYS_biomass',-0.000309143);
model = setParam(model,'eq','pMET_biomass',-6.37596E-05);
model = setParam(model,'eq','pCYS_biomass',-0.00014299);
model = setParam(model,'eq','pPHE_biomass',-0.000273721);
model = setParam(model,'eq','pTYR_biomass',-0.000125466);
model = setParam(model,'eq','pTHR_biomass',-0.000243847);
model = setParam(model,'eq','pTRP_biomass',-5.79283E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000273125);
model = setParam(model,'eq','pARG_biomass',-0.000304441);
model = setParam(model,'eq','pALA_biomass',-0.000371609);
model = setParam(model,'eq','pASP_biomass',-0.000378715);
model = setParam(model,'eq','pASN_biomass',-0.000381983);
model = setParam(model,'eq','pGLU_biomass',-0.000448159);
model = setParam(model,'eq','pGLN_biomass',-0.000514553);
model = setParam(model,'eq','pGLY_biomass',-0.000412796);
model = setParam(model,'eq','pPRO_biomass',-0.000530388);
model = setParam(model,'eq','pSER_biomass',-0.000381706);

save([data '/mat/wr1_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.000279294);
model = setParam(model,'eq','sGLC_biomass',-1.14182E-05);
model = setParam(model,'eq','sFRU_biomass',1.81288E-05);
model = setParam(model,'eq','sSUCROSE_biomass',-8.28742E-05);
model = setParam(model,'eq','Raffinose_biomass',-4.98775E-05);
model = setParam(model,'eq','Stachyose_biomass',-9.42738E-05);
model = setParam(model,'eq','Verbascose_biomass',0);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-2.22688E-07);
model = setParam(model,'eq','Palmitate_biomass',-0.000257001);
model = setParam(model,'eq','Stearate_biomass',-9.58998E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000399619);
model = setParam(model,'eq','Linoleate_biomass',-0.000794742);
model = setParam(model,'eq','Linolenate_biomass',-0.000140808);
model = setParam(model,'eq','Arachidate_biomass',-8.16673E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000565487);
model = setParam(model,'eq','CELLWALL_biomass',-0.004892237);
model = setParam(model,'eq','pHIS_biomass',-0.000120535);
model = setParam(model,'eq','pILE_biomass',-0.000289592);
model = setParam(model,'eq','pLEU_biomass',-0.000568932);
model = setParam(model,'eq','pLYS_biomass',-0.000407265);
model = setParam(model,'eq','pMET_biomass',-8.39967E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000188374);
model = setParam(model,'eq','pPHE_biomass',-0.0003606);
model = setParam(model,'eq','pTYR_biomass',-0.000165289);
model = setParam(model,'eq','pTHR_biomass',-0.000321244);
model = setParam(model,'eq','pTRP_biomass',-7.63146E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000359815);
model = setParam(model,'eq','pARG_biomass',-0.00040107);
model = setParam(model,'eq','pALA_biomass',-0.000489557);
model = setParam(model,'eq','pASP_biomass',-0.000498917);
model = setParam(model,'eq','pASN_biomass',-0.000503224);
model = setParam(model,'eq','pGLU_biomass',-0.000590403);
model = setParam(model,'eq','pGLN_biomass',-0.000677871);
model = setParam(model,'eq','pGLY_biomass',-0.000543816);
model = setParam(model,'eq','pPRO_biomass',-0.000698732);
model = setParam(model,'eq','pSER_biomass',-0.000502858);

save([data '/mat/wr2_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
model = setParam(model,'eq','Starch_biomass',0.001688528);
model = setParam(model,'eq','sGLC_biomass',9.70032E-05);
model = setParam(model,'eq','sFRU_biomass',0.000152151);
model = setParam(model,'eq','sSUCROSE_biomass',0.000138938);
model = setParam(model,'eq','Raffinose_biomass',-0.000115188);
model = setParam(model,'eq','Stachyose_biomass',-0.000273148);
model = setParam(model,'eq','Verbascose_biomass',-1.2172E-06);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-6.99735E-08);
model = setParam(model,'eq','Palmitate_biomass',-0.000182233);
model = setParam(model,'eq','Stearate_biomass',-5.21744E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000333979);
model = setParam(model,'eq','Linoleate_biomass',-0.000771494);
model = setParam(model,'eq','Linolenate_biomass',-0.000153908);
model = setParam(model,'eq','Arachidate_biomass',-4.75287E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000499537);
model = setParam(model,'eq','CELLWALL_biomass',-0.000533425);
model = setParam(model,'eq','pHIS_biomass',-4.48256E-05);
model = setParam(model,'eq','pILE_biomass',-0.000107696);
model = setParam(model,'eq','pLEU_biomass',-0.000211579);
model = setParam(model,'eq','pLYS_biomass',-0.000151457);
model = setParam(model,'eq','pMET_biomass',-3.12374E-05);
model = setParam(model,'eq','pCYS_biomass',-7.0054E-05);
model = setParam(model,'eq','pPHE_biomass',-0.000134103);
model = setParam(model,'eq','pTYR_biomass',-6.14689E-05);
model = setParam(model,'eq','pTHR_biomass',-0.000119467);
model = setParam(model,'eq','pTRP_biomass',-2.83805E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000133811);
model = setParam(model,'eq','pARG_biomass',-0.000149153);
model = setParam(model,'eq','pALA_biomass',-0.00018206);
model = setParam(model,'eq','pASP_biomass',-0.000185541);
model = setParam(model,'eq','pASN_biomass',-0.000187143);
model = setParam(model,'eq','pGLU_biomass',-0.000219564);
model = setParam(model,'eq','pGLN_biomass',-0.000252092);
model = setParam(model,'eq','pGLY_biomass',-0.000202238);
model = setParam(model,'eq','pPRO_biomass',-0.00025985);
model = setParam(model,'eq','pSER_biomass',-0.000187007);

save([data '/mat/wr3_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',0.001401472);
model = setParam(model,'eq','sGLC_biomass',5.71886E-05);
model = setParam(model,'eq','sFRU_biomass',5.92674E-05);
model = setParam(model,'eq','sSUCROSE_biomass',0.000365063);
model = setParam(model,'eq','Raffinose_biomass',-6.85502E-05);
model = setParam(model,'eq','Stachyose_biomass',-0.000201478);
model = setParam(model,'eq','Verbascose_biomass',-2.20024E-06);
model = setParam(model,'eq','Maltose_biomass',0);
model = setParam(model,'eq','Myristate_biomass',-1.64076E-08);
model = setParam(model,'eq','Palmitate_biomass',-9.56995E-05);
model = setParam(model,'eq','Stearate_biomass',-1.48818E-05);
model = setParam(model,'eq','Oleate_biomass',-0.000138924);
model = setParam(model,'eq','Linoleate_biomass',-0.00050047);
model = setParam(model,'eq','Linolenate_biomass',-8.979E-05);
model = setParam(model,'eq','Arachidate_biomass',-1.21132E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000280331);
model = setParam(model,'eq','CELLWALL_biomass',0);
model = setParam(model,'eq','pHIS_biomass',-2.21235E-05);
model = setParam(model,'eq','pILE_biomass',-5.31528E-05);
model = setParam(model,'eq','pLEU_biomass',-0.000104424);
model = setParam(model,'eq','pLYS_biomass',-7.47509E-05);
model = setParam(model,'eq','pMET_biomass',-1.54171E-05);
model = setParam(model,'eq','pCYS_biomass',-3.45749E-05);
model = setParam(model,'eq','pPHE_biomass',-6.61859E-05);
model = setParam(model,'eq','pTYR_biomass',-3.03377E-05);
model = setParam(model,'eq','pTHR_biomass',-5.89623E-05);
model = setParam(model,'eq','pTRP_biomass',-1.40071E-05);
model = setParam(model,'eq','pVAL_biomass',-6.60417E-05);
model = setParam(model,'eq','pARG_biomass',-7.36139E-05);
model = setParam(model,'eq','pALA_biomass',-8.98552E-05);
model = setParam(model,'eq','pASP_biomass',-9.15732E-05);
model = setParam(model,'eq','pASN_biomass',-9.23636E-05);
model = setParam(model,'eq','pGLU_biomass',-0.000108365);
model = setParam(model,'eq','pGLN_biomass',-0.000124419);
model = setParam(model,'eq','pGLY_biomass',-9.9814E-05);
model = setParam(model,'eq','pPRO_biomass',-0.000128248);
model = setParam(model,'eq','pSER_biomass',-9.22966E-05);

save([data '/mat/wr4_c2-c3.mat'],'model')


