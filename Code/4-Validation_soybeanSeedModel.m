%Setting different biomass composition according to the conditions each
%plant was submitted in the field:

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

% Setting Raven solver
% getpref('RAVEN','solver')
% setRavenSolver('gurobi')


%% Load curated model 'Glycine max' GmaxGEM with RAVEN toolbox
% The curation was obtained with manualCuration.m script

% Load the curated model
load([data '/mat/soybeanSeedModel.mat'],'model')

% Firsly, we set all biomass composition to either 0 or fixed as
% experimental data provided. First setting all biomass to zero:
BiomassRxns = model.rxns(endsWith(model.rxns,'_biomass'));
model = setParam(model, 'eq', BiomassRxns,0);

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then, setting the different biomass composition.
% Define the model constraints from Allen and Young, 2013:
model = setParam(model, 'eq', 'Starch_biomass', 0.0);
model = setParam(model, 'eq', 'sGLC_biomass', 0.0);
model = setParam(model, 'eq', 'sFRU_biomass', 0.0);
model = setParam(model, 'eq', 'sSUCROSE_biomass', 0.0);
model = setParam(model, 'eq', 'Raffinose_biomass', 0.0);
model = setParam(model, 'eq', 'Stachyose_biomass', 0.0);
model = setParam(model, 'eq', 'Verbascose_biomass', 0.0);
model = setParam(model, 'eq', 'Myristate_biomass', 0.0);
model = setParam(model, 'eq', 'Palmitate_biomass', 0.0);
model = setParam(model, 'eq', 'Stearate_biomass', 0.0);
model = setParam(model, 'eq', 'Oleate_biomass', 0.0);
model = setParam(model, 'eq', 'Linoleate_biomass', 0.0);
model = setParam(model, 'eq', 'Linolenate_biomass', 0.0);
model = setParam(model, 'eq', 'Arachidate_biomass', 0.0);
model = setParam(model, 'eq', 'Cellulose_biomass', -3.2564);
model = setParam(model, 'eq', 'pHIS_biomass', -0.1092);
model = setParam(model, 'eq', 'pILE_biomass', -0.1634);
model = setParam(model, 'eq', 'pLEU_biomass', -0.2488);
model = setParam(model, 'eq', 'pLYS_biomass', -0.1339);
model = setParam(model, 'eq', 'pMET_biomass', -0.0765);
model = setParam(model, 'eq', 'pCYS_biomass', -0.0685);
model = setParam(model, 'eq', 'pPHE_biomass', -0.2202);
model = setParam(model, 'eq', 'pTYR_biomass', -0.02);
model = setParam(model, 'eq', 'pTHR_biomass', -0.1146);
model = setParam(model, 'eq', 'pTRP_biomass', -0.0);
model = setParam(model, 'eq', 'pVAL_biomass', -0.2088);
model = setParam(model, 'eq', 'pARG_biomass', -0.3098);
model = setParam(model, 'eq', 'pALA_biomass', -0.1908);
model = setParam(model, 'eq', 'pASP_biomass', -0.1274);
model = setParam(model, 'eq', 'pASN_biomass', -0.1274);
model = setParam(model, 'eq', 'pGLU_biomass', -0.15775);
model = setParam(model, 'eq', 'pGLN_biomass', -0.15775);
model = setParam(model, 'eq', 'pGLY_biomass', -0.3325);
model = setParam(model, 'eq', 'pPRO_biomass', -0.184);
model = setParam(model, 'eq', 'pSER_biomass', -0.1016);
model = setParam(model, 'eq', 'sMAL_biomass', 0.0);
model = setParam(model, 'eq', 'sCIT_biomass', 0.0);
model = setParam(model, 'eq', 'sASP_biomass', 0.0);
model = setParam(model, 'eq', 'sALA_biomass', 0.0);
model = setParam(model, 'eq', 'sASN_biomass', 0.0);
model = setParam(model, 'eq', 'sGLU_biomass', 0.0);
model = setParam(model, 'eq', 'sGLN_biomass', 0.0);
model = setParam(model, 'eq', 'sGLY_biomass', 0.0);
model = setParam(model, 'eq', 'sTRP_biomass', 0.0);
model = setParam(model, 'eq', 'sTYR_biomass', 0.0);
model = setParam(model, 'eq', 'sPHE_biomass', 0.0);
model = setParam(model, 'eq', 'sVAL_biomass', 0.0);
model = setParam(model, 'eq', 'sILE_biomass', 0.0);
model = setParam(model, 'eq', 'sLEU_biomass', 0.0);
model = setParam(model, 'eq', 'sMET_biomass', 0.0);
model = setParam(model, 'eq', 'sCYS_biomass', 0.0);
model = setParam(model, 'eq', 'sHIS_biomass', 0.0);
model = setParam(model, 'eq', 'sLYS_biomass', 0.0);
model = setParam(model, 'eq', 'sISOCIT_biomass', 0.0);
model = setParam(model, 'eq', 'sFUM_biomass', 0.0);
model = setParam(model, 'eq', 'sUREA_biomass', 0.0);
model = setParam(model, 'eq', 'sGLYCOLATE_biomass', 0.0);
model = setParam(model, 'eq', 'sSER_biomass', 0.0);
model = setParam(model, 'eq', 'sTHR_biomass', 0.0);
model = setParam(model, 'eq', 'ASCORBATE_biomass', 0.0);
model = setParam(model, 'eq', 's2KG_biomass', 0.0);
model = setParam(model, 'eq', 'sPYROGLU_biomass', 0.0);
model = setParam(model, 'eq', 'NH4_tx', 0.0);
model = setParam(model, 'eq', 'ATPase_tx', 0.0);
model = setParam(model, 'eq', 'sSUC_biomass', 0.0);
model = setParam(model, 'eq', 'Galacturonate_biomass', 0.0);
model = setParam(model, 'eq', 'Maltose_biomass', 0.0);
model = setParam(model, 'eq', 'PAPS_biomass', 0.0);
model = setParam(model, 'eq', 'CO_A_biomass', 0.0);
model = setParam(model, 'eq', 'GALACTOSE_biomass', 0.0);
model = setParam(model, 'eq', 'Glycerol_biomass', -0.15);
model = setParam(model, 'eq', 'ARABINOSE_biomass', 0.0);
model = setParam(model, 'eq', 'Xylan_biomass', 0.0);
model = setParam(model, 'eq', 'TRANSCINNAMATE_biomass', 0.0);
model = setParam(model, 'eq', 'sPRO_biomass', 0.0);
model = setParam(model, 'eq', 'RHAMNOSE_biomass', 0.0);
model = setParam(model, 'eq', '2PHENYLETHANOL_biomass', 0.0);
model = setParam(model, 'eq', 'UMP_biomass', 0.0);
model = setParam(model, 'eq', 'CMP_biomass', 0.0);
model = setParam(model, 'eq', 'AMP_biomass', 0.0);
model = setParam(model, 'eq', 'GMP_biomass', 0.0);
model = setParam(model, 'eq', 'CoumOL_biomass', 0.0);
model = setParam(model, 'eq', 'INDOLE_ACETATE_AUXIN_biomass', 0.0);
model = setParam(model, 'eq', 'sETOH_biomass', 0.0);
model = setParam(model, 'eq', 'Glycerol_3P_biomass', 0.0);
model = setParam(model, 'eq', 'CHLOROPHYLL_A_biomass', 0.0);
model = setParam(model, 'eq', 'CHLOROPHYLL_B_biomass', 0.0);
model = setParam(model, 'eq', 'ConOL_biomass', 0.0);
model = setParam(model, 'eq', 'sGABA_biomass', 0.0);
model = setParam(model, 'eq', 'Dodecanoate_biomass', 0.0);
model = setParam(model, 'eq', 'dTMP_biomass', 0.0);
model = setParam(model, 'eq', 'dCMP_biomass', 0.0);
model = setParam(model, 'eq', 'dAMP_biomass', 0.0);
model = setParam(model, 'eq', 'SinapOL_biomass', 0.0);
model = setParam(model, 'eq', 'GLUTATHIONE_biomass', 0.0);
model = setParam(model, 'eq', 'TAG_biomass', 0.0);
model = setParam(model, 'eq', 'FUCOSE_biomass', 0.0);
model = setParam(model, 'eq', 'FattyAcid_biomass', -0.4925);
model = setParam(model, 'eq', 'FAD_biomass', 0.0);
model = setParam(model, 'eq', 'THF_biomass', 0.0);
model = setParam(model, 'eq', 'sORN_biomass', 0.0);
model = setParam(model, 'eq', 'SAM_biomass', 0.0);
model = setParam(model, 'eq', 'dGMP_biomass', 0.0);
model = setParam(model, 'eq', 'sARG_biomass', 0.0);
model = setParam(model, 'eq', 'Methanetiol_biomass', 0.0);
model = setParam(model, 'eq', 'NADPHoxc_tx', 0.0);
model = setParam(model, 'eq', 'NADPHoxm_tx', 0.0);
model = setParam(model, 'eq', 'NADPHoxp_tx', 0.0);
model = setParam(model, 'eq', 'ASN_tx', 0.3297);
model = setParam(model, 'eq', 'GLN_tx', 1.6475);
model = setParam(model, 'eq', 'GLC_tx', 1.4153);
model = setParam(model, 'eq', 'Sucrose_tx', 2.8335);


% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx';
          'Photon_tx';
          'NO3_tx';
          'SO4_tx'};

model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

save([data '/mat/allen2013.mat'],'model')
exportToExcelFormat(model, [root '/scrap/soybeanSeedModel_validation.xlsx']);
% load([data '/mat/allen2013.mat'],'model')

% Set one reaction as objective function
model = setParam(model,'obj','Photon_tx',-1);
          
% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y);
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/validation_minPhot.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/validation__minPhot_exch.txt']);

load([data '/mat/allen2013.mat'],'model')

% Set one reaction as objective function
model = setParam(model,'obj','ATPase_tx',1);
          
% Check objective
x = find(model.c==1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y);
% Check constraints
% printConstraints(model)

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/validation_maxATP.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/validation__maxATP_exch.txt']);

% checking specific flux:
ListOfRxnFlxs = {'SO4_tx';'NO3_tx';'GLC_tx';'H2O_tx';'O2_tx';'Sucrose_tx';
    'Photon_tx';'PROTON_tx';'CO2_tx';'ATPase_tx';'Ca_tx';'Fe_tx';'K_tx';
    'Mg_tx';'NADPHoxc_tx';'NADPHoxm_tx';'NADPHoxp_tx';'NADPHoxx_tx';'NH4_tx';'Starch_tx';
    'TAG_tx';'Pi_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



