%% PROTOCOL FOR THE FBA WITH SOYBEAN 'Glycine max' GENOME-SCALE METABOLIC MODEL (SoyGEM) USING THE RAVEN TOOLBOX 

% This protocol was prepared to run FBA with 'Glycine max', available on GitHub

% Define paths for model reconstruction.
% clear; % clean workspace
clc; % clean command window
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

% Setting Raven solver
% setRavenSolver('gurobi')
% getpref('RAVEN','solver')


%% Running pFBA with the soymodel specific to seed(embryo)model
% With script Biomass_soybeanSeedModel.m we obtained 8 models that differed 
% between then according to treatment (wr1, wr2, wr3 and wr4) or period 
% of the harvest (the first:c1-c2 or the second: c2-c3).
% Then, separetely, each model was used to apply pFBA analysis with RAVEN toolbox.
% Firsly,  all biomass components were set either to 0 or fixed as
% experimental data provided.
% Then, we run pFBA for every different biomass composition
% Additionally, the Photon influx was minimized as part of the objective function

% ########################################################################################################
% The light input was set in orther to reach the carbon conversion efficiencies of 85%
% as experimentally calculated in previous embryo studies (Allen et al., 2009, 2013)
% Thus it was differently set for each treatment, as follow:
% 70-82DAP
%  WR1	Photon_tx = 0.0096
%  WR2	Photon_tx = 0.0346
%  WR3	Photon_tx = 0.0480
%  WR4	Photon_tx = 0.0590
% 
% 82-97DAP
%  WR1	Photon_tx = 0.0636
%  WR2	Photon_tx = 0.1248
%  WR3	Photon_tx = 0.1576
%  WR4	Photon_tx = 0.1256

%% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
load([data '/mat/wr1_c1-c2.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.0096);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.2806E-04);
model = setParam(model, 'eq','NADPHoxc_tx',-5.8673E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-5.8673E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-5.8673E-05);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
              'ASN_tx'; % <=> ASN_c[c]	asparagine
              'GLN_tx'; % <=> L-glutamine[c]
              'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
              'NH4_tx'; % <=> H+[c] + ammonia[c]
                 }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_00096_WR1_C1-C2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_00096_WR1_C1-C2_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    % 'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    % 'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    % 'NADPHoxp_tx';'NADPHoxx_tx'
    };
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};

ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% Running sampling 
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-7.7928e-07,10);
model = setParam(model, 'var','CO2_tx',-0.0041786,10);
model = setParam(model, 'var','Glycerol_biomass',-7.5491e-05,10);
model = setParam(model, 'var','H2O_tx',-0.0063777,10);
model = setParam(model, 'var','Linoleate_biomass',-9.5098e-05,10);
model = setParam(model, 'var','Linolenate_biomass',-3.9044e-05,10);
model = setParam(model, 'var','Myristate_biomass',-6.0377e-08,10);
model = setParam(model, 'var','O2_tx',0.001945,10);
model = setParam(model, 'var','Oleate_biomass',-4.1153e-05,10);
model = setParam(model, 'var','PROTON_tx',0.0038257,10);
model = setParam(model, 'var','Palmitate_biomass',-3.9609e-05,10);
model = setParam(model, 'var','Photon_tx',0.0096,10);
model = setParam(model, 'var','Raffinose_biomass',-3.9432e-06,10);
model = setParam(model, 'var','Stachyose_biomass',-5.5356e-06,10);
model = setParam(model, 'var','Starch_biomass',-0.0004561,10);
model = setParam(model, 'var','Stearate_biomass',-9.4174e-06,10);
model = setParam(model, 'var','sFRU_biomass',-4.1257e-05,10);
model = setParam(model, 'var','sGLC_biomass',-2.5783e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',-0.00014212,10);
model = setParam(model, 'var','GLN_tx',0.0010959,10);
model = setParam(model, 'var','CELLWALL_biomass',-0.0012329,10);

% Running Sampling with RAVEN
goodRxns = [];
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_00096_WR1_C1-C2_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_00096_WR1_C1-C2_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_00096_WR1_C1-C2_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_00096_WR1_C1-C2_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_00096_WR1_C1-C2_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_00096_WR1_C1-C2_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_00096_WR1_C1-C2_NADPH_median.xls');



%% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr2_c1-c2.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.0346);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.3901E-04);
model = setParam(model, 'eq','NADPHoxc_tx',-5.9890E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-5.9890E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-5.9890E-05);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_00346_WR2_C1-C2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_00346_WR2_C1-C2_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    % 'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    % 'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    % 'NADPHoxp_tx';'NADPHoxx_tx'
    };
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% Running sampling 
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-3.2881e-06,10);
model = setParam(model, 'var','CO2_tx',-0.012066,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00028972,10);
model = setParam(model, 'var','H2O_tx',-0.017819,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00035836,10);
model = setParam(model, 'var','Linolenate_biomass',-0.00010296,10);
model = setParam(model, 'var','Myristate_biomass',-1.711e-07,10);
model = setParam(model, 'var','O2_tx',0.0040709,10);
model = setParam(model, 'var','Oleate_biomass',-0.00020708,10);
model = setParam(model, 'var','PROTON_tx',0.010206,10);
model = setParam(model, 'var','Palmitate_biomass',-0.00014947,10);
model = setParam(model, 'var','Photon_tx',0.0346,10);
model = setParam(model, 'var','Raffinose_biomass',-2.5423e-06,10);
model = setParam(model, 'var','Starch_biomass',-0.0012812,10);
model = setParam(model, 'var','Stearate_biomass',-4.2623e-05,10);
model = setParam(model, 'var','sFRU_biomass',-9.6664e-05,10);
model = setParam(model, 'var','sGLC_biomass',-6.1786e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',-0.00044332,10);
model = setParam(model, 'var','GLN_tx',0.0030373,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.0032362,10);

% Running Sampling with RAVEN
goodRxns = [];
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_00346_WR2_C1-C2_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_00346_WR2_C1-C2_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_00346_WR2_C1-C2_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_00346_WR2_C1-C2_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_00346_WR2_C1-C2_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_00346_WR2_C1-C2_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_00346_WR2_C1-C2_NADPH_median.xls');

%% #########################################################################################################
%  >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr3_c1-c2.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.0480);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.4544E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-2.7271E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-2.7271E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-2.7271E-04);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_00480_WR3_C1-C2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_00480_WR3_C1-C2_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% Running sampling 
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-5.719e-06,10);
model = setParam(model, 'var','CO2_tx',-0.021082,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00039641,10);
model = setParam(model, 'var','H2O_tx',-0.031871,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00043201,10);
model = setParam(model, 'var','Linolenate_biomass',-8.8513e-05,10);
model = setParam(model, 'var','Myristate_biomass',-3.0714e-07,10);
model = setParam(model, 'var','O2_tx',0.0094419,10);
model = setParam(model, 'var','Oleate_biomass',-0.00036444,10);
model = setParam(model, 'var','PROTON_tx',0.019111,10);
model = setParam(model, 'var','Palmitate_biomass',-0.00022167,10);
model = setParam(model, 'var','Photon_tx',0.0480,10);
model = setParam(model, 'var','Raffinose_biomass',-1.1359e-05,10);
model = setParam(model, 'var','Starch_biomass',-0.0024534,10);
model = setParam(model, 'var','Stearate_biomass',-6.8068e-05,10);
model = setParam(model, 'var','sFRU_biomass',-0.00015148,10);
model = setParam(model, 'var','sGLC_biomass',-0.00011493,10);
model = setParam(model, 'var','sSUCROSE_biomass',-0.00073554,10);
model = setParam(model, 'var','GLN_tx',0.0054016,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.0064541,10);

% Running Sampling with RAVEN
goodRxns = [];
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_00480_WR3_C1-C2_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_00480_WR3_C1-C2_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_00480_WR3_C1-C2_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_00480_WR3_C1-C2_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_00480_WR3_C1-C2_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_00480_WR3_C1-C2_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_00480_WR3_C1-C2_NADPH_median.xls');

%% #########################################################################################################
%  >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr4_c1-c2.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.0590);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.0100E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-1.1222E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-1.1222E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-1.1222E-04);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

      
% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_00590_WR4_C1-C2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_00590_WR4_C1-C2_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% Running sampling

% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-7.2845e-06,10);
model = setParam(model, 'var','CO2_tx',-0.022333,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00048237,10);
model = setParam(model, 'var','H2O_tx',-0.033699,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00053685,10);
model = setParam(model, 'var','Linolenate_biomass',-9.6743e-05,10);
model = setParam(model, 'var','Myristate_biomass',-3.2293e-07,10);
model = setParam(model, 'var','O2_tx',0.0085128,10);
model = setParam(model, 'var','Oleate_biomass',-0.00045974,10);
model = setParam(model, 'var','PROTON_tx',0.020574,10);
model = setParam(model, 'var','Palmitate_biomass',-0.00025269,10);
model = setParam(model, 'var','Photon_tx',0.0590,10);
model = setParam(model, 'var','Raffinose_biomass',-1.0234e-05,10);
model = setParam(model, 'var','Starch_biomass',-0.0020904,10);
model = setParam(model, 'var','Stearate_biomass',-8.2698e-05,10);
model = setParam(model, 'var','sFRU_biomass',-0.00012202,10);
model = setParam(model, 'var','sGLC_biomass',-0.00012883,10);
model = setParam(model, 'var','sSUCROSE_biomass',-0.00069067,10);
model = setParam(model, 'var','GLN_tx',0.0058966,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.0068832,10);

% Running Sampling with RAVEN
goodRxns = []; 
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_00590_WR4_C1-C2_sampling.csv');


% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_00590_WR4_C1-C2_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_00590_WR4_C1-C2_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_00590_WR4_C1-C2_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_00590_WR4_C1-C2_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_00590_WR4_C1-C2_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_00590_WR4_C1-C2_NADPH_median.xls');

%% #########################################################################################################
%  >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr1_c2-c3.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.0636);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.9163E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-6.5737E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-6.5737E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-6.5737E-04);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_00636_WR1_C2-C3.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_00636_WR1_C2-C3_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ListOfRxnFlxs = {'GLUTAMATE_DEHYDROGENASE_NADP_RXN_c';'H2O_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% Running sampling
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-5.1193e-06,10);
model = setParam(model, 'var','CO2_tx',-0.014203,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00039955,10);
model = setParam(model, 'var','H2O_tx',-0.020866,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00051034,10);
model = setParam(model, 'var','Linolenate_biomass',-0.00010418,10);
model = setParam(model, 'var','Myristate_biomass',-1.9923e-07,10);
model = setParam(model, 'var','O2_tx',0.0034385,10);
model = setParam(model, 'var','Oleate_biomass',-0.00033052,10);
model = setParam(model, 'var','PROTON_tx',0.011953,10);
model = setParam(model, 'var','Palmitate_biomass',-0.00016865,10);
model = setParam(model, 'var','Photon_tx',0.0636,10);
model = setParam(model, 'var','Raffinose_biomass',3.3073e-06,10);
model = setParam(model, 'var','Stachyose_biomass',4.4285e-06,10);
model = setParam(model, 'var','Starch_biomass',-0.0012895,10);
model = setParam(model, 'var','Stearate_biomass',-7.2506e-05,10);
model = setParam(model, 'var','sFRU_biomass',-1.0749e-05,10);
model = setParam(model, 'var','sGLC_biomass',-3.6789e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',-0.00047488,10);
model = setParam(model, 'var','GLN_tx',0.0036148,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.0038158,10);

% Running Sampling with RAVEN
goodRxns = []; 
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;

writetable(struct2table(out),'photon_00636_WR1_C2-C3_sampling.csv');


% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_00636_WR1_C2-C3_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_00636_WR1_C2-C3_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_00636_WR1_C2-C3_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_00636_WR1_C2-C3_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_00636_WR1_C2-C3_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_00636_WR1_C2-C3_NADPH_median.xls');


%% #########################################################################################################
%  >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr2_c2-c3.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.1248);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.5809E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-1.7565E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-1.7565E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-1.7565E-03);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_01248_WR2_C2-C3.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_01248_WR2_C2-C3_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% Running sampling 
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-8.1667e-06,10);
model = setParam(model, 'var','CO2_tx',-0.016528,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00056949,10);
model = setParam(model, 'var','H2O_tx',-0.024518,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00079474,10);
model = setParam(model, 'var','Linolenate_biomass',-0.00014081,10);
model = setParam(model, 'var','Myristate_biomass',-2.2269e-07,10);
model = setParam(model, 'var','O2_tx',0.0013926,10);
model = setParam(model, 'var','Oleate_biomass',-0.00039962,10);
model = setParam(model, 'var','PROTON_tx',0.015044,10);
model = setParam(model, 'var','Palmitate_biomass',-0.000257,10);
model = setParam(model, 'var','Photon_tx',0.1248,10);
model = setParam(model, 'var','Raffinose_biomass',-4.9877e-05,10);
model = setParam(model, 'var','Stachyose_biomass',-9.4274e-05,10);
model = setParam(model, 'var','Starch_biomass',-0.00031033,10);
model = setParam(model, 'var','Stearate_biomass',-9.59e-05,10);
model = setParam(model, 'var','sFRU_biomass',1.8129e-05,10);
model = setParam(model, 'var','sGLC_biomass',1.1418e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',8.2874e-05,10);
model = setParam(model, 'var','GLN_tx',0.0047621,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.0043602,10);

% Running Sampling with RAVEN
goodRxns = [];
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_01248_WR2_C2-C3_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_01248_WR2_C2-C3_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_01248_WR2_C2-C3_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_01248_WR2_C2-C3_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_01248_WR2_C2-C3_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_01248_WR2_C2-C3_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_01248_WR2_C2-C3_NADPH_median.xls');


%% #########################################################################################################
%  >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

load([data '/mat/wr3_c2-c3.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.1576);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9864E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3182E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3182E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3182E-03);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_01576_WR3_C2-C3.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_01576_WR3_C2-C3_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% Running sampling 
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-4.7529e-06,10);
model = setParam(model, 'var','CO2_tx',-0.0071198,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00050286,10);
model = setParam(model, 'var','H2O_tx',-0.0085866,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00077149,10);
model = setParam(model, 'var','Linolenate_biomass',-0.00015391,10);
model = setParam(model, 'var','Myristate_biomass',-6.9974e-08,10);
model = setParam(model, 'var','O2_tx',-0.0048815,10);
model = setParam(model, 'var','Oleate_biomass',-0.00033398,10);
model = setParam(model, 'var','PROTON_tx',0.0038902,10);
model = setParam(model, 'var','Palmitate_biomass',-0.00018223,10);
model = setParam(model, 'var','Photon_tx',0.1576,10);
model = setParam(model, 'var','Raffinose_biomass',-0.00011519,10);
model = setParam(model, 'var','Stachyose_biomass',-0.00027315,10);
model = setParam(model, 'var','Starch_biomass',0.0018761,10);
model = setParam(model, 'var','Stearate_biomass',-5.2174e-05,10);
model = setParam(model, 'var','Verbascose_biomass',-1.2172e-06,10);
model = setParam(model, 'var','sFRU_biomass',0.00015215,10);
model = setParam(model, 'var','sGLC_biomass',9.7003e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',0.00013894,10);
model = setParam(model, 'var','GLN_tx',0.001771,10);
model = setParam(model, 'var', 'CELLWALL_biomass',-0.00065308,10);

% Running Sampling with RAVEN
goodRxns = [];
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;
writetable(struct2table(out),'photon_01576_WR3_C2-C3_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_01576_WR3_C2-C3_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_01576_WR3_C2-C3_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_01576_WR3_C2-C3_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_01576_WR3_C2-C3_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_01576_WR3_C2-C3_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_01576_WR3_C2-C3_NADPH_median.xls');

%% #########################################################################################################
%   >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
load([data '/mat/wr4_c2-c3.mat'],'model')

% Closing all exchange rxns (_tx)
exchangeRxns = model.rxns(endsWith(model.rxns,'_tx'));
model = setParam(model, 'eq', exchangeRxns, 0);

% Then only a set of reactions will be openned:
organicSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
                'ASN_tx'; % <=> ASN_c[c]	asparagine
                'GLN_tx' % <=> L-glutamine[c]
                  }; 
model = setParam(model, 'lb', organicSource,0);
model = setParam(model, 'ub', organicSource,1000);

% Photon influx:
model = setParam(model, 'eq', 'Photon_tx',0.1256);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9829E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3144E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3144E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3144E-03);

% Carbon dioxide, oxygen, water, and proton were allowed freely exchange with the environment
freeXc = {'O2_tx';  %  <=> oxygen[c]
          'CO2_tx';  %  <=> CO2[c]
          'H2O_tx'; %  <=> H2O[c]
          'PROTON_tx'}; %  <=> H+[c]
model = setParam(model,'lb',freeXc,-1000);
model = setParam(model,'ub',freeXc,1000);

% nitrate was set as inorganic nitrogen source.
nitroSource = {'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
               'NH4_tx'; % <=> H+[c] + ammonia[c]
               };  
model = setParam(model,'lb',nitroSource,0);
model = setParam(model,'ub',nitroSource,1000);

% Mineral nutrients reactions were set to be opened
mineraNuts = {'SO4_tx';  % 3 Pumped-PROTON_c[c] => 3 H+[c] + sulfate[c]
               'Ca_tx';  % Pumped-PROTON_c[c] => H+[c] + Ca+2_c[c]
               'K_tx';  %  Pumped-PROTON_c[c] => H+[c] + K+_c[c]
               'Pi_tx';  %   3 Pumped-PROTON_c[c] => 3 H+[c] + phosphate[c]
               'Mg_tx';  %   Pumped-PROTON_c[c] => H+[c] + Mg+2_c[c]
               'TRANS-RXN0-200_EXP_1'; %  => H+[c] + Zn2+[c]
               'Fe_tx'}; % <=> IRON(II)_c[c]
model = setParam(model,'lb', mineraNuts, 0);
model = setParam(model,'ub', mineraNuts, 1000);

% Set one reaction as objective function
minSource = {'Sucrose_tx'; % Pumped-PROTON_c[c] => H+[c] + sucrose[c]
             'ASN_tx'; % <=> ASN_c[c]	asparagine
             'GLN_tx'; % <=> L-glutamine[c]
             'NO3_tx';  %   2 Pumped-PROTON_c[c] => 2 H+[c] + nitrate[c]
             'NH4_tx'; % <=> H+[c] + ammonia[c]
                }; 
model = setParam(model,'obj',minSource,-1);

% Check objective
x = find(model.c==-1); % choose ==1 or ==-1
y = model.rxns(x);
disp(y)
% Check constraints
% printConstraints(model)

%  Check model integrity
% checkModelStruct(model, true);

% Check fluxes
sol=solveLP(model,1);
printFluxes(model,sol.x, true)
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_01256_WR4_C2-C3.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_01256_WR4_C2-C3_toSampl.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ListOfRxnFlxs = {'Ca_tx';'Cellulose_biomass';'ConOL_biomass';'CoumOL_biomass';'ALA_tx'};
ListOfRxnFlxs = {'GLYCEROL_KIN_RXN_c';'Mitochondrial_ATP_Synthase_m';'PEPDEPHOS_RXN_c';'ACETATE__COA_LIGASE_RXN_p';'ACETATE__COA_LIGASE_RXN_x';'ACETYL_COA_CARBOXYLTRANSFER_RXN_c'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% Running sampling
% Adding a variation =10% of the exchanges rxns within each model
model = setParam(model, 'var','Arachidate_biomass',-1.2113e-06,10);
model = setParam(model, 'var','CO2_tx',-0.0029852,10);
model = setParam(model, 'var','Glycerol_biomass',-0.00028221,10);
model = setParam(model, 'var','H2O_tx',-0.0033626,10);
model = setParam(model, 'var','Linoleate_biomass',-0.00050047,10);
model = setParam(model, 'var','Linolenate_biomass',-8.979e-05,10);
model = setParam(model, 'var','Myristate_biomass',-1.6408e-08,10);
model = setParam(model, 'var','O2_tx',-0.0036618,10);
model = setParam(model, 'var','Oleate_biomass',-0.00013892,10);
model = setParam(model, 'var','PROTON_tx',0.0015401,10);
model = setParam(model, 'var','Palmitate_biomass',-9.57e-05,10);
model = setParam(model, 'var','Photon_tx',0.1256,10);
model = setParam(model, 'var','Raffinose_biomass',-6.855e-05,10);
model = setParam(model, 'var','Stachyose_biomass',-0.00020148,10);
model = setParam(model, 'var','Starch_biomass',0.0015572,10);
model = setParam(model, 'var','Stearate_biomass',-1.4882e-05,10);
model = setParam(model, 'var','Verbascose_biomass',-2.2002e-06,10);
model = setParam(model, 'var','sFRU_biomass',5.9267e-05,10);
model = setParam(model, 'var','sGLC_biomass',5.7189e-05,10);
model = setParam(model, 'var','sSUCROSE_biomass',0.00036506,10);
model = setParam(model, 'var','GLN_tx',0.00087405,10);

% Running Sampling with RAVEN
goodRxns = []; 
sol = randomSampling(model,10000,true,true,true,goodRxns,true);
sol = full(sol);

% output: ec
out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;

writetable(struct2table(out),'photon_01256_WR4_C2-C3_sampling.csv');

% getting the MEAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outATP_mean.rxns    = model.rxns(rxnIdx);
    outATP_mean.rxnNames= model.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end
writecell(outATP_mean,'photon_01256_WR4_C2-c3_ATP_mean.xls');

clear x y outATP_mean rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADH_mean.rxns    = model.rxns(rxnIdx);
    outNADH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end
writecell(outNADH_mean,'photon_01256_WR4_C2-c3_NADH_mean.xls');

clear x y outNADH_mean rxnIdx 
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.mean,true);
    outNADPH_mean.rxns    = model.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end
writecell(outNADPH_mean,'photon_01256_WR4_C2-c3_NADPH_mean.xls');

% getting the MEDIAN of the fluxes of all reactions with ATP, NADH (NAD) and NADPH (NADP):
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outATP_median.rxns    = model.rxns(rxnIdx);
    outATP_median.rxnNames= model.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end
writecell(outATP_median,'photon_01256_WR4_C2-C3_ATP_median.xls');

clear x y outATP_median rxnIdx 
for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADH_median.rxns    = model.rxns(rxnIdx);
    outNADH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end
writecell(outNADH_median,'photon_01256_WR4_C2-C3_NADH_median.xls');

clear x y outNADH_median rxnIdx
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,out.median,true);
    outNADPH_median.rxns    = model.rxns(rxnIdx);
    outNADPH_median.rxnNames= model.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end
writecell(outNADPH_median,'photon_01256_WR4_C2-C3_NADPH_median.xls');


% save([data '/mat/alles_FBA_soybean_mitPhoton.mat'])

