%% Manual curation
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

%% Load 'Glycine max' GmaxGEM (*.yml) with RAVEN toolbox
% Source: Moreira et al., 2019. Plant Physiology, Vol. 180, pp. 1912–1929.
% Some corrections were made with correctSoybean-GEM.m and with
% replace_genesID2.m

% Load the correct model
load([data '/mat/SoyGEM_correctGeneID.mat'],'model');
% or open it with readYAMLmodel
% model = readYAMLmodel([data '/templateModels/SoyGEM_correctGeneID.yml']);
 
% ----------------------------------------------
%% Removing duplicated and incorrected rxns
% 1st: import in 'HOME' the list of rxns to be removed as a cell array: 'RxnsToRmvSoy.txt'
% then remove all of them:
size(model.rxns) % = 3001
modeltmp = removeReactions(model,RxnsToRmvSoy,true, true, true);
size(modeltmp.rxns) % = 2947

%% keeping only one of each duplicated rxns and only the corrected one:
% 2nd: adding rxns-> one of each duplicated and deleted before; missing rxns;
% and uncorrected rxns:
fid = fopen('RxnsToAddSoy.txt');
loadedData = textscan(fid,'%q %q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);
clear rxnsToAdd
rxnsToAdd.equations     = regexprep(loadedData{1},'***','');
rxnsToAdd.rxnNames      = regexprep(loadedData{2},'***','');
rxnsToAdd.grRules       = regexprep(loadedData{3},'***','');
rxnsToAdd.subSystems    = regexprep(loadedData{4},'***','');
rxnsToAdd.eccodes       = regexprep(loadedData{5},'***','');
rxnsToAdd.rxns          = rxnsToAdd.rxnNames;

size(modeltmp.rxns) % = 2947
modeltmp2 = addRxns(modeltmp,rxnsToAdd,3,'',true,true);
size(modeltmp2.rxns) % = 2984
model = modeltmp2;

%% Some metabolites/rxns had no name, these now have the same name as their identifier.
emptyNames=cellfun(@isempty,model.rxnNames);
model.rxnNames(emptyNames)=model.rxns(emptyNames);

% Changing specific bounds
model = setParam(model, 'ub', 'HYDROXYPYRUVATE_REDUCTASE_RXN_NADP_c',0);
model = setParam(model, 'lb', 'HYDROXYPYRUVATE_REDUCTASE_RXN_NADP_c',-1000);

% Standardize NADPH and NADP name from plastid and mitochondria
% Replace NADPH_p with NADPH, and replace NADPH_m with NADPH
% model.mets = regexprep(model.mets,'NADPH_p','NADPH'); 
% model.mets = regexprep(model.mets,'NADPH_m','NADPH'); 
% model.metNames = regexprep(model.metNames,'NADPH_p','NADPH');
% model.metNames = regexprep(model.metNames,'NADPH_m','NADPH');
% model.grRules = regexprep(model.grRules,'NADPH_p','NADPH'); 
% model.grRules = regexprep(model.grRules,'NADPH_m','NADPH');
% 
% model.mets = regexprep(model.mets,'NADP_p','NADP+'); 
% model.mets = regexprep(model.mets,'NADP_m','NADP+'); 
% model.mets = regexprep(model.mets,'NAD_c','NAD+'); 
% model.metNames = regexprep(model.metNames,'NADP_p','NADP+');
% model.metNames = regexprep(model.metNames,'NADP_m','NADP+');
% model.metNames = regexprep(model.metNames,'NAD_c','NAD+');
% model.grRules = regexprep(model.grRules,'NADP_p','NADP+'); 
% model.grRules = regexprep(model.grRules,'NADP_m','NADP+'); 
% model.grRules = regexprep(model.grRules,'NAD_c','NAD+'); 


% check and saving model
save([data '/mat/soybeanSeedModel.mat'],'model')
exportToExcelFormat(model, [root '/scrap/soybeanSeedModel.xlsx']);
writeYAMLmodel(model, [data '/templateModels/soybeanSeedModel.yml']);
exportModel(model,[data '/templateModels/soybeanSeedModel.xml']);

