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
model = setParam(model,'eq','Starch_biomass',-0.000456096);
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
model = setParam(model,'eq','Glycerol_biomass',-7.54909E-05);
model = setParam(model,'eq','CELLWALL_biomass',-0.00123291);
model = setParam(model,'eq','pHIS_biomass',-2.45217E-05);
model = setParam(model,'eq','pILE_biomass',-5.7501E-05);
model = setParam(model,'eq','pLEU_biomass',-0.000112966);
model = setParam(model,'eq','pLYS_biomass',-8.21869E-05);
model = setParam(model,'eq','pMET_biomass',-1.6999E-05);
model = setParam(model,'eq','pCYS_biomass',-3.69116E-05);
model = setParam(model,'eq','pPHE_biomass',-7.39453E-05);
model = setParam(model,'eq','pTYR_biomass',-3.42605E-05);
model = setParam(model,'eq','pTHR_biomass',-6.27591E-05);
model = setParam(model,'eq','pTRP_biomass',-1.6015E-05);
model = setParam(model,'eq','pVAL_biomass',-7.00837E-05);
model = setParam(model,'eq','pARG_biomass',-8.27645E-05);
model = setParam(model,'eq','pALA_biomass',-8.99037E-05);
model = setParam(model,'eq','pASP_biomass',-9.92929E-05);
model = setParam(model,'eq','pASN_biomass',-0.000100033);
model = setParam(model,'eq','pGLU_biomass',-0.000135875);
model = setParam(model,'eq','pGLN_biomass',-0.00013679);
model = setParam(model,'eq','pGLY_biomass',-9.51432E-05);
model = setParam(model,'eq','pPRO_biomass',-0.000135665);
model = setParam(model,'eq','pSER_biomass',-9.5906E-05);

save([data '/mat/wr1_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.001281192);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000289717);
model = setParam(model,'eq','CELLWALL_biomass',-0.003236166);
model = setParam(model,'eq','pHIS_biomass',-6.79593E-05);
model = setParam(model,'eq','pILE_biomass',-0.000159358);
model = setParam(model,'eq','pLEU_biomass',-0.000313075);
model = setParam(model,'eq','pLYS_biomass',-0.000227772);
model = setParam(model,'eq','pMET_biomass',-4.71108E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000102297);
model = setParam(model,'eq','pPHE_biomass',-0.000204931);
model = setParam(model,'eq','pTYR_biomass',-9.49492E-05);
model = setParam(model,'eq','pTHR_biomass',-0.00017393);
model = setParam(model,'eq','pTRP_biomass',-4.43839E-05);
model = setParam(model,'eq','pVAL_biomass',-0.00019423);
model = setParam(model,'eq','pARG_biomass',-0.000229373);
model = setParam(model,'eq','pALA_biomass',-0.000249158);
model = setParam(model,'eq','pASP_biomass',-0.00027518);
model = setParam(model,'eq','pASN_biomass',-0.000277231);
model = setParam(model,'eq','pGLU_biomass',-0.000376563);
model = setParam(model,'eq','pGLN_biomass',-0.0003791);
model = setParam(model,'eq','pGLY_biomass',-0.000263679);
model = setParam(model,'eq','pPRO_biomass',-0.00037598);
model = setParam(model,'eq','pSER_biomass',-0.000265793);

save([data '/mat/wr2_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.002453398);
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
model = setParam(model,'eq','Linoleate_biomass',-0.00043201 );
model = setParam(model,'eq','Linolenate_biomass',-8.85133E-05);
model = setParam(model,'eq','Arachidate_biomass',-5.71899E-06);
model = setParam(model,'eq','Glycerol_biomass',-0.000396412);
model = setParam(model,'eq','CELLWALL_biomass',-0.00645414);
model = setParam(model,'eq','pHIS_biomass',-0.000120862);
model = setParam(model,'eq','pILE_biomass',-0.000283409);
model = setParam(model,'eq','pLEU_biomass',-0.000556785);
model = setParam(model,'eq','pLYS_biomass',-0.00040508);
model = setParam(model,'eq','pMET_biomass',-8.37839E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000181929);
model = setParam(model,'eq','pPHE_biomass',-0.000364459);
model = setParam(model,'eq','pTYR_biomass',-0.000168862);
model = setParam(model,'eq','pTHR_biomass',-0.000309325);
model = setParam(model,'eq','pTRP_biomass',-7.89343E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000345426);
model = setParam(model,'eq','pARG_biomass',-0.000407927);
model = setParam(model,'eq','pALA_biomass',-0.000443114);
model = setParam(model,'eq','pASP_biomass',-0.000489391);
model = setParam(model,'eq','pASN_biomass',-0.00049304 );
model = setParam(model,'eq','pGLU_biomass',-0.000669696);
model = setParam(model,'eq','pGLN_biomass',-0.000674207);
model = setParam(model,'eq','pGLY_biomass',-0.000468938);
model = setParam(model,'eq','pPRO_biomass',-0.000668658);
model = setParam(model,'eq','pSER_biomass',-0.000472698);

save([data '/mat/wr3_c1-c2.mat'],'model')

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.002090365);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000482367);
model = setParam(model,'eq','CELLWALL_biomass',-0.00688320);
model = setParam(model,'eq','pHIS_biomass',-0.000131935);
model = setParam(model,'eq','pILE_biomass',-0.000309375);
model = setParam(model,'eq','pLEU_biomass',-0.000607799);
model = setParam(model,'eq','pLYS_biomass',-0.000442194);
model = setParam(model,'eq','pMET_biomass',-9.14604E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000198597);
model = setParam(model,'eq','pPHE_biomass',-0.000397851);
model = setParam(model,'eq','pTYR_biomass',-0.000184333);
model = setParam(model,'eq','pTHR_biomass',-0.000337666);
model = setParam(model,'eq','pTRP_biomass',-8.61664E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000377075);
model = setParam(model,'eq','pARG_biomass',-0.000445302);
model = setParam(model,'eq','pALA_biomass',-0.000483713);
model = setParam(model,'eq','pASP_biomass',-0.00053423 );
model = setParam(model,'eq','pASN_biomass',-0.000538213);
model = setParam(model,'eq','pGLU_biomass',-0.000731055);
model = setParam(model,'eq','pGLN_biomass',-0.000735979);
model = setParam(model,'eq','pGLY_biomass',-0.000511904);
model = setParam(model,'eq','pPRO_biomass',-0.000729922);
model = setParam(model,'eq','pSER_biomass',-0.000516008);

save([data '/mat/wr4_c1-c2.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.001289543);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000399553);
model = setParam(model,'eq','CELLWALL_biomass',-0.003815785);
model = setParam(model,'eq','pHIS_biomass',-8.08805E-05);
model = setParam(model,'eq','pILE_biomass',-0.000189657);
model = setParam(model,'eq','pLEU_biomass',-0.0003726);
model = setParam(model,'eq','pLYS_biomass',-0.000271079);
model = setParam(model,'eq','pMET_biomass',-5.6068E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000121746);
model = setParam(model,'eq','pPHE_biomass',-0.000243895);
model = setParam(model,'eq','pTYR_biomass',-0.000113002);
model = setParam(model,'eq','pTHR_biomass',-0.000207);
model = setParam(model,'eq','pTRP_biomass',-5.28227E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000231159);
model = setParam(model,'eq','pARG_biomass',-0.000272984);
model = setParam(model,'eq','pALA_biomass',-0.000296531);
model = setParam(model,'eq','pASP_biomass',-0.0003275);
model = setParam(model,'eq','pASN_biomass',-0.000329941);
model = setParam(model,'eq','pGLU_biomass',-0.000448159);
model = setParam(model,'eq','pGLN_biomass',-0.000451178);
model = setParam(model,'eq','pGLY_biomass',-0.000313813);
model = setParam(model,'eq','pPRO_biomass',-0.000447465);
model = setParam(model,'eq','pSER_biomass',-0.000316328);

save([data '/mat/wr1_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',-0.000310327);
model = setParam(model,'eq','sGLC_biomass',1.14182E-05);
model = setParam(model,'eq','sFRU_biomass',1.81288E-05);
model = setParam(model,'eq','sSUCROSE_biomass',8.28742E-05);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000569493);
model = setParam(model,'eq','CELLWALL_biomass',-0.004360169);
model = setParam(model,'eq','pHIS_biomass',-0.000106552);
model = setParam(model,'eq','pILE_biomass',-0.000249853);
model = setParam(model,'eq','pLEU_biomass',-0.000490861);
model = setParam(model,'eq','pLYS_biomass',-0.000357118);
model = setParam(model,'eq','pMET_biomass',-7.38639E-05);
model = setParam(model,'eq','pCYS_biomass',-0.000160388);
model = setParam(model,'eq','pPHE_biomass',-0.000321307);
model = setParam(model,'eq','pTYR_biomass',-0.000148868);
model = setParam(model,'eq','pTHR_biomass',-0.000272701);
model = setParam(model,'eq','pTRP_biomass',-6.95884E-05);
model = setParam(model,'eq','pVAL_biomass',-0.000304528);
model = setParam(model,'eq','pARG_biomass',-0.000359628);
model = setParam(model,'eq','pALA_biomass',-0.000390649);
model = setParam(model,'eq','pASP_biomass',-0.000431447);
model = setParam(model,'eq','pASN_biomass',-0.000434664);
model = setParam(model,'eq','pGLU_biomass',-0.000590403);
model = setParam(model,'eq','pGLN_biomass',-0.00059438);
model = setParam(model,'eq','pGLY_biomass',-0.000413416);
model = setParam(model,'eq','pPRO_biomass',-0.000589489);
model = setParam(model,'eq','pSER_biomass',-0.00041673);

save([data '/mat/wr2_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
model = setParam(model,'eq','Starch_biomass',0.001876142);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000502861);
model = setParam(model,'eq','CELLWALL_biomass',-0.00065308);
model = setParam(model,'eq','pHIS_biomass',-3.96253E-05);
model = setParam(model,'eq','pILE_biomass',-9.29172E-05);
model = setParam(model,'eq','pLEU_biomass',-0.000182545);
model = setParam(model,'eq','pLYS_biomass',-0.000132808);
model = setParam(model,'eq','pMET_biomass',-2.74691E-05);
model = setParam(model,'eq','pCYS_biomass',-5.96464E-05);
model = setParam(model,'eq','pPHE_biomass',-0.00011949 );
model = setParam(model,'eq','pTYR_biomass',-5.53623E-05);
model = setParam(model,'eq','pTHR_biomass',-0.000101414);
model = setParam(model,'eq','pTRP_biomass',-2.58791E-05);
model = setParam(model,'eq','pVAL_biomass',-0.00011325);
model = setParam(model,'eq','pARG_biomass',-0.000133741);
model = setParam(model,'eq','pALA_biomass',-0.000145278);
model = setParam(model,'eq','pASP_biomass',-0.00016045);
model = setParam(model,'eq','pASN_biomass',-0.000161646);
model = setParam(model,'eq','pGLU_biomass',-0.000219564);
model = setParam(model,'eq','pGLN_biomass',-0.000221043);
model = setParam(model,'eq','pGLY_biomass',-0.000153744);
model = setParam(model,'eq','pPRO_biomass',-0.000219224);
model = setParam(model,'eq','pSER_biomass',-0.000154977);

save([data '/mat/wr3_c2-c3.mat'],'model')

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
model = setParam(model, 'eq', BiomassRxns,0);

model = setParam(model,'eq','Starch_biomass',0.00155719);
model = setParam(model,'eq','sGLC_biomass',5.71886E-05);
model = setParam(model,'eq','sFRU_biomass',5.92674E-05);
model = setParam(model,'eq','sSUCROSE_biomass',0.00036506);
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
model = setParam(model,'eq','Glycerol_biomass',-0.000282207);
model = setParam(model,'eq','CELLWALL_biomass',0); % 0.00005045);
model = setParam(model,'eq','pHIS_biomass',-1.95569E-05);
model = setParam(model,'eq','pILE_biomass',-4.5859E-05);
model = setParam(model,'eq','pLEU_biomass',-9.00946E-05);
model = setParam(model,'eq','pLYS_biomass',-6.55469E-05);
model = setParam(model,'eq','pMET_biomass',-1.35573E-05);
model = setParam(model,'eq','pCYS_biomass',-2.94382E-05);
model = setParam(model,'eq','pPHE_biomass',-5.89739E-05);
model = setParam(model,'eq','pTYR_biomass',-2.73239E-05);
model = setParam(model,'eq','pTHR_biomass',-5.00525E-05);
model = setParam(model,'eq','pTRP_biomass',-1.27725E-05);
model = setParam(model,'eq','pVAL_biomass',-5.58942E-05);
model = setParam(model,'eq','pARG_biomass',-6.60075E-05);
model = setParam(model,'eq','pALA_biomass',-7.17012E-05);
model = setParam(model,'eq','pASP_biomass',-7.91894E-05);
model = setParam(model,'eq','pASN_biomass',-7.97798E-05);
model = setParam(model,'eq','pGLU_biomass',-0.000108365);
model = setParam(model,'eq','pGLN_biomass',-0.000109095);
model = setParam(model,'eq','pGLY_biomass',-7.58799E-05);
model = setParam(model,'eq','pPRO_biomass',-0.000108197);
model = setParam(model,'eq','pSER_biomass',-7.64882E-05);

save([data '/mat/wr4_c2-c3.mat'],'model')


