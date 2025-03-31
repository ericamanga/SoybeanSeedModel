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
%  WR1	Photon_tx = 0.010350  
%  WR2	Photon_tx = 0.036335
%  WR3	Photon_tx = 0.051372
%  WR4	Photon_tx = 0.062182
% 
% 82-97DAP
%  WR1	Photon_tx = 0.064948
%  WR2	Photon_tx = 0.122706
%  WR3	Photon_tx = 0.141885
%  WR4	Photon_tx = 0.112092

%% #############################################################################################################################
%% ###########################################  MAINTENANCE X10   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.2806E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-5.8673E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-5.8673E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-5.8673E-04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_MaintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_MainX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.3901E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-5.9890E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-5.9890E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-5.9890E-04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_maintX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.4544E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-2.7271E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-2.7271E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-2.7271E-03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_maintX10.txt']);


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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.0100E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-1.1222E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-1.1222E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-1.1222E-03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_maintX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.9163E-02);
model = setParam(model, 'eq','NADPHoxc_tx',-6.5737E-03);
model = setParam(model, 'eq','NADPHoxm_tx',-6.5737E-03);
model = setParam(model, 'eq','NADPHoxp_tx',-6.5737E-03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_maintX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.5809E-01);
model = setParam(model, 'eq','NADPHoxc_tx',-1.7565E-02);
model = setParam(model, 'eq','NADPHoxm_tx',-1.7565E-02);
model = setParam(model, 'eq','NADPHoxp_tx',-1.7565E-02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_maintX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9864E-01);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3182E-02);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3182E-02);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3182E-02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_maintX10.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9829E-01);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3144E-02);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3144E-02);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3144E-02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_maintX10.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_maintX10.txt']);

% #########################################################################################################
%% ###########################################  MAINTENANCE X0.1   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.2806E-05);
model = setParam(model, 'eq','NADPHoxc_tx',-5.8673E-06);
model = setParam(model, 'eq','NADPHoxm_tx',-5.8673E-06);
model = setParam(model, 'eq','NADPHoxp_tx',-5.8673E-06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_MaintX01.txt']); % 01 refers to 0.1


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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.3901E-05);
model = setParam(model, 'eq','NADPHoxc_tx',-5.9890E-06);
model = setParam(model, 'eq','NADPHoxm_tx',-5.9890E-06);
model = setParam(model, 'eq','NADPHoxp_tx',-5.9890E-06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_MaintX01.txt']); % 01 refers to 0.1


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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.4544E-04);
model = setParam(model, 'eq','NADPHoxc_tx',-2.7271E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-2.7271E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-2.7271E-05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_MaintX01.txt']); % 01 refers to 0.1

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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.0100E-04);
model = setParam(model, 'eq','NADPHoxc_tx',-1.1222E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-1.1222E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-1.1222E-05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_MaintX01.txt']); % 01 refers to 0.1



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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-5.9163E-04);
model = setParam(model, 'eq','NADPHoxc_tx',-6.5737E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-6.5737E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-6.5737E-05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_MaintX01.txt']); % 01 refers to 0.1

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-1.5809E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-1.7565E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-1.7565E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-1.7565E-04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_MaintX01.txt']); % 01 refers to 0.1

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9864E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3182E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3182E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3182E-04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_MaintX01.txt']); % 01 refers to 0.1



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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-2.9829E-03);
model = setParam(model, 'eq','NADPHoxc_tx',-3.3144E-04);
model = setParam(model, 'eq','NADPHoxm_tx',-3.3144E-04);
model = setParam(model, 'eq','NADPHoxp_tx',-3.3144E-04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_MaintX01.txt']); % 01 refers to 0.1
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_MaintX01.txt']); % 01 refers to 0.1

% #########################################################################################################
%% ###########################################  MAINTENANCE X2   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00105612);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000117346);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000117346);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000117346);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00107802);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00011978);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00011978);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00011978);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0049088);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00054542);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00054542);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00054542);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_maintX2.txt']);


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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00202);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00022444);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00022444);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00022444);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0118326);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00131474);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00131474);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00131474);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.031618);
model = setParam(model, 'eq','NADPHoxc_tx',-0.003513);
model = setParam(model, 'eq','NADPHoxm_tx',-0.003513);
model = setParam(model, 'eq','NADPHoxp_tx',-0.003513);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.059728);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0066364);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0066364);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0066364);


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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_maintX2.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.059658);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0066288);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0066288);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0066288);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_maintX2.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_maintX2.txt']);

% #########################################################################################################
%% ###########################################  MAINTENANCE X0.5   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00026403);
model = setParam(model, 'eq','NADPHoxc_tx',-2.93365E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-2.93365E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-2.93365E-05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.000269505);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000029945);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000029945);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000029945);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0012272);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000136355);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000136355);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000136355);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_maintX05.txt']);


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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.000505);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00005611);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00005611);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00005611);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00295815);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000328685);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000328685);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000328685);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0079045);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00087825);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00087825);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00087825);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.014932);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0016591);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0016591);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0016591);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_maintX05.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0149145);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0016572);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0016572);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0016572);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_maintX05.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_maintX05.txt']);

% #########################################################################################################
%% ###########################################  MAINTENANCE X5   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0026403);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000293365);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000293365);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000293365);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00269505);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00029945);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00029945);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00029945);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.012272);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00136355);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00136355);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00136355);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_maintX5.txt']);


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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00505);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0005611);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0005611);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0005611);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0295815);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00328685);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00328685);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00328685);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.079045);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0087825);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0087825);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0087825);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.14932);
model = setParam(model, 'eq','NADPHoxc_tx',-0.016591);
model = setParam(model, 'eq','NADPHoxm_tx',-0.016591);
model = setParam(model, 'eq','NADPHoxp_tx',-0.016591);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_maintX5.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.149145);
model = setParam(model, 'eq','NADPHoxc_tx',-0.016572);
model = setParam(model, 'eq','NADPHoxm_tx',-0.016572);
model = setParam(model, 'eq','NADPHoxp_tx',-0.016572);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_maintX5.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_maintX5.txt']);

% #########################################################################################################
%% ###########################################  MAINTENANCE X0.2   ##############################################################
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
model = setParam(model, 'eq', 'Photon_tx',0.010350);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.000105612);
model = setParam(model, 'eq','NADPHoxc_tx',-1.17346E-05);
model = setParam(model, 'eq','NADPHoxm_tx',-1.17346E-05);
model = setParam(model, 'eq','NADPHoxp_tx',-1.17346E-05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0010350_WR1_C1-C2_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.036335);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.000107802);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000011978);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000011978);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000011978);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0036335_WR2_C1-C2_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.051372);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00049088);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000054542);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000054542);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000054542);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0051372_WR3_C1-C2_toSampl_maintX02.txt']);


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
model = setParam(model, 'eq', 'Photon_tx',0.062182);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.000202);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000022444);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000022444);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000022444);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0062182_WR4_C1-C2_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.064948);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.00118326);
model = setParam(model, 'eq','NADPHoxc_tx',-0.000131474);
model = setParam(model, 'eq','NADPHoxm_tx',-0.000131474);
model = setParam(model, 'eq','NADPHoxp_tx',-0.000131474);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0064948_WR1_C2-C3_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.122706);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0031618);
model = setParam(model, 'eq','NADPHoxc_tx',-0.0003513);
model = setParam(model, 'eq','NADPHoxm_tx',-0.0003513);
model = setParam(model, 'eq','NADPHoxp_tx',-0.0003513);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0122706_WR2_C2-C3_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.141885);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0059728);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00066364);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00066364);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00066364);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0141885_WR3_C2-C3_toSampl_maintX02.txt']);

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
model = setParam(model, 'eq', 'Photon_tx',0.112092);

% Fixed ATP and NAPH maintenance (values from pair of cotiledones at 
% 24-48h in Moreira et al., 2019)
model = setParam(model, 'eq','ATPase_tx',-0.0059658);
model = setParam(model, 'eq','NADPHoxc_tx',-0.00066288);
model = setParam(model, 'eq','NADPHoxm_tx',-0.00066288);
model = setParam(model, 'eq','NADPHoxp_tx',-0.00066288);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_maintX02.txt']);
printFluxes(model, sol.x,true,'',[root '/outputs/FBA_photon_0112092_WR4_C2-C3_toSampl_maintX02.txt']);

