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

%% ########################################################################################################
% The light input was set in orther to find the RUBISCO activation
% Thus it was differently set for each treatment, as follow:
% C1-C2 and C2-C3 (x2)
% WR1, WR2, WR3 and WR4 (x4)
% with photon_tx ranging from 0.01 to 0.15 (15 photon inputs x 8 treatments= 120 simlations)

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.01 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
load([data '/mat/wr1_c1-c2.mat'],'model')
disp(model)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.01);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_001_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.02 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR3_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR4_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.02);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_002_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.03 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.03);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_003_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.04 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.04);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_004_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.05 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.05);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_005_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.06 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.06);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_006_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.07 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.07);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_007_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.08 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.08);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_008_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.09 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.09);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_009_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.10 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR1_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR2_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR3_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR4_C1-C2.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR1_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR3_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.10);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_010_WR4_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.11 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.11);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_011_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.12 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.12);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_012_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 
% ########################################################################################################
% ------------------------------------------ PHOTON = 0.13 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.13);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_013_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.14 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.14);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_014_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% ########################################################################################################
% ------------------------------------------ PHOTON = 0.15 --------------------------------------------------- 
% #########################################################################################################
% >>>>> For ΔC1-C2 in treatment WR 1 (x 1000)
clear; % clean workspace
cd ../;
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean.m file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR1_C1-C2.txt']);

% % checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

 
%% >>>>> For ΔC1-C2 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR2_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC1-C2 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR3_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC1-C2 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR4_C1-C2.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])
% 

 
%% >>>>> For ΔC2-C3 in treatment WR 1 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR1_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])

%% >>>>> For ΔC2-C3 in treatment WR 2 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR2_C2-C3.txt']);

% checking specific flux:
ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
    'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
    'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
    'NADPHoxp_tx';'NADPHoxx_tx'};
ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
sol.x(ResultsFlx)
cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])



%% >>>>> For ΔC2-C3 in treatment WR 3 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR3_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])


%% >>>>> For ΔC2-C3 in treatment WR 4 (x 1000)
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
model = setParam(model, 'eq', 'Photon_tx',0.15);

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
printFluxes(model, sol.x,false,'',[root '/outputs/FBA_photon_015_WR4_C2-C3.txt']);

% checking specific flux:
% ListOfRxnFlxs = {'Photon_tx';'Sucrose_tx';'ASN_tx';'GLN_tx';%'GLC_tx';'FRU_tx';'GLU_tx';'ALA_tx';
%     'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p';'O2_tx';'CO2_tx';...
%     'H2O_tx';'NO3_tx'; 'NH4_tx';'ATPase_tx';'NADPHoxc_tx';'NADPHoxm_tx'; ...
%     'NADPHoxp_tx';'NADPHoxx_tx'};
% ResultsFlx = getIndexes(model,ListOfRxnFlxs,'rxns');
% sol.x(ResultsFlx)
% cell2table([num2cell(sol.x(ResultsFlx)), ListOfRxnFlxs])





