%% PROTOCOL FOR THE CORRECTIONS OF THE MODEL SOYBEAN 'Glycine max'
% GENOME-SCALE METABOLIC MODEL (SoyGEM) USING THE RAVEN TOOLBOX 
% Source: Moreira et al., 2019. Plant Physiology, Vol. 180, pp. 1912–1929.

% Define paths for model reconstruction.
% clear
clc; % clean command window
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];

%% Correcting it to run into RAVEN
modelSoy = importModel([data '/templateModels/File1_soyModel.xml'],false,false,true);

modelSoy.id = 'gmx';
modelSoy.name = 'gmx';

%% Fix empty compartment annotation
% Make compartments
modelSoy.comps={'c';'m';'p';'v';'x'};
modelSoy.compNames={'cytosol';'mitochondrion';...
                    'plastid';'vacuole';...
                    'peroxisome'};
% Assign metabolites to the correct compartment
suffix=endsWith(modelSoy.mets,'_c');
modelSoy.metComps(suffix) = 1;
suffix=endsWith(modelSoy.mets,'_m');
modelSoy.metComps(suffix) = 2;
suffix=endsWith(modelSoy.mets,'_p');
modelSoy.metComps(suffix) = 3;
suffix=endsWith(modelSoy.mets,'_v');
modelSoy.metComps(suffix) = 4;
suffix=endsWith(modelSoy.mets,'_x');
modelSoy.metComps(suffix) = 5;

checkModelStruct(modelSoy,false)

% Metabolites with no name: fill with mets field, remove compartment ID
hasNoMetName = find(cellfun(@isempty,modelSoy.metNames));
newMetNames = modelSoy.mets(hasNoMetName);
newMetNames = regexprep(newMetNames,'_(.)_','');
modelSoy.metNames(hasNoMetName) = newMetNames;

% Check model again
checkModelStruct(modelSoy,false)

% Replace duplicate metabolites
% All metabolites named arachidoyl-CoA, by position, metabolite ID and compartment
pos = find(strcmp('arachidoyl-CoA',modelSoy.metNames));
metID = modelSoy.mets(strcmp('arachidoyl-CoA',modelSoy.metNames));
comp = modelSoy.comps(modelSoy.metComps(pos));
[num2cell(pos), metID, comp]

% Include coefficients from 1019 with 669
modelSoy.S(669,:) = modelSoy.S(669,:) + modelSoy.S(1019,:);
modelSoy = removeMets(modelSoy,'CPD-9965_c',false,true,true,true);

% All metabolites named esculin, by position, metabolite ID and compartment
pos = find(strcmp('esculin',modelSoy.metNames));
metID = modelSoy.mets(strcmp('esculin',modelSoy.metNames));
comp = modelSoy.comps(modelSoy.metComps(pos));
[num2cell(pos), metID, comp]
% CPD-11682 does not exist on MetaCyc
% Include coefficients from 1317 with 1704
modelSoy.S(1704,:) = modelSoy.S(1704,:) + modelSoy.S(1317,:);
modelSoy = removeMets(modelSoy,'CPD-11682_c',false,true,true,true);

% All metabolites named L-tyrosyl-tRNA<SUP>tyr</SUP>, by position, metabolite ID and compartment
pos = find(strcmp('L-tyrosyl-tRNA<SUP>tyr</SUP>',modelSoy.metNames));
metID = modelSoy.mets(strcmp('L-tyrosyl-tRNA<SUP>tyr</SUP>',modelSoy.metNames));
comp = modelSoy.comps(modelSoy.metComps(pos));
[num2cell(pos), metID, comp]
% Wrongly named, should be valine instead of tyrosine
modelSoy.metNames{2809} = 'L-valyl-tRNA<SUP>val</SUP>';

% Correct chemical formula (molybdenum is Mo, not MO)
pos = getIndexes(modelSoy,'CPD-3_c','mets');
modelSoy.metFormulas(pos) = {'Mo1O4'};

%% Reaction, metabolite and gene identifiers all have characters that are
% not allowed in identifiers in an SBML file. These characters include
% anything that is not a letter, number or underscore, for instance +, -
% and .. When exporting the model, they are now replaced by place-holders
% (e.g. _DASH_ for -), to avoid such illegal characters ending up in the
% SBML file. For simplicity (and faster model load/export), we can remove
% and/or replace some of these characters:

% This is not essential, if you keep these special characters in the MATLAB
% structure they will just be replaced during each round of loading and
% exporting the model.

% Replace - and . with _: (note that we use \. in the regular expression
% instead of ., as . is otherwise a special character in regexp).
modelSoy.rxns = regexprep(modelSoy.rxns,'-|\.','_'); % Replace - and . with _
modelSoy.mets = regexprep(modelSoy.mets,'-|\.','_'); % Replace - and . with _
modelSoy.genes = regexprep(modelSoy.genes,'-|\.','_'); % Replace - and . with _

% Make sure that any changes applied to .genes field are also applied to
% the .grRules field, otherwise they do not match anymore.
modelSoy.grRules = regexprep(modelSoy.grRules,'-|\.','_'); % Replace - and . with _
modelSoy.rxns = regexprep(modelSoy.rxns,'+',''); % Remove +
modelSoy.mets = regexprep(modelSoy.mets,'+',''); % Remove +
modelSoy.rxns = regexprep(modelSoy.rxns,'-',''); % Remove -
modelSoy.mets = regexprep(modelSoy.mets,'-',''); % Remove -

%% additional clean of the model
modelSoy.metNames=regexprep(modelSoy.metNames,'<sub>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</sub>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'<SUB>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</SUB>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'<SUP>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</SUP>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'<sup>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</sup>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'<i>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</i>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'<I>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'</I>','');
modelSoy.metNames=regexprep(modelSoy.metNames,'&alpha','alpha');
modelSoy.metNames=regexprep(modelSoy.metNames,'&Alpha','alpha');
modelSoy.metNames=regexprep(modelSoy.metNames,'&beta','beta');
modelSoy.metNames=regexprep(modelSoy.metNames,'&gamma','gamma');
modelSoy.metNames=regexprep(modelSoy.metNames,'&delta','delta');
modelSoy.metNames=regexprep(modelSoy.metNames,'&Delta','delta');
modelSoy.metNames=regexprep(modelSoy.metNames,'&zeta','zeta');
modelSoy.metNames=regexprep(modelSoy.metNames,'&omega','omega');
modelSoy.metNames=regexprep(modelSoy.metNames,'&rarr','rarr');
modelSoy.metNames=regexprep(modelSoy.metNames,'&epsilon','epsilon');
modelSoy.metNames=regexprep(modelSoy.metNames,'&prime','prime');
modelSoy.metNames=regexprep(modelSoy.metNames,'&lambda','lambda');

modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&mdash','-mdash');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&alpha','alpha');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&Alpha','alpha');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&beta','beta');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&gamma','gamma');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&delta','delta');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&Delta','delta');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&zeta','zeta');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&omega','omega');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&rarr','rarr');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&epsilon','epsilon');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&prime','prime');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'&lambda','lambda');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<em>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</em>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<i>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</i>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<I>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</I>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<small>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</small>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<sub>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</sub>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<SUB>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</SUB>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<SUP>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</SUP>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'<sup>','');
modelSoy.rxnNames=regexprep(modelSoy.rxnNames,'</sup>','');

checkModelStruct(modelSoy,false)

% Set arbitrary reaction as objective function
modelSoy = setParam(modelSoy,'obj','Photon_tx',1);

% Saving the model
exportModel(modelSoy,[data '/templateModels/SoyGEM_correct.xml']);
exportToExcelFormat(modelSoy, [root '/scrap/SoyGEM_correct.xlsx']);
save([data '/mat/SoyGEM_correct.mat'],'modelSoy')

