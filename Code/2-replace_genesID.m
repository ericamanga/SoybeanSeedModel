% Script to replace genes ID from a csv file to the soybean model

clear; % clean workspace
clc; % clean command window
if ~exist([pwd() '/FBA_soybean.m']); error(['Make sure that '...
        'your Current Folder is the one containing the FBA_soybean file.']); end
cd ../;
root = [pwd() '/SoybeanSeedModel'];
data = [root '/data/'];
code = [root '/code/'];
cd(data)

%% Replacing the genes ID in the model with informations from csv file
% Loading the metabolic model from an XML file after correctSoybean-GEM.m
model = importModel([data '/templateModels/SoyGEM_correct.xml']);

model.id = 'soy';
disp(model)
 
% Make use of COBRA functions generateRules and creategrRulesField, to
% avoid having to replacing all gene names in grRules. First, make rules
% field based on original gene names.
model = generateRules(model);

% Load the gene ID replacement information from a CSV file
geneReplacementData = readtable([data '/templateModels/soy_gene_mapping.csv']);

% Check for the correct column names in the CSV file
if ~ismember('OriginalGeneID', geneReplacementData.Properties.VariableNames) || ...
   ~ismember('ReplacementGeneID', geneReplacementData.Properties.VariableNames)
    error("CSV file does not contain the required column names 'OriginalGeneID' and 'ReplacementGeneID'");
end

% Loop through the geneShortNames field and replace gene IDs, give genes the
% same value
for i = 1:length(model.geneShortNames)
    [~, idx] = ismember(model.geneShortNames{i}, geneReplacementData.OriginalGeneID);
    if idx > 0
        model.genes{i} = geneReplacementData.ReplacementGeneID{idx};
        model.geneShortNames{i} = geneReplacementData.ReplacementGeneID{idx};
    end
end

% Show the first 10 entries of each field
model.genes(1:10)
model.geneShortNames(1:10)

% As we have just renamed the genes and not changed their position in the
% genes field, we can rebuild a new grRules field, based on the previously
% defined rules field, and now using the new gene names.
model = creategrRulesField(model);
model = rmfield(model,'rules'); % Not essential to do this, but cleanest to just remove the COBRA field.
% Show the first 10 entries
model.grRules(1:10)

% Saving and exporting the updated model
save([data '/mat/SoyGEM_correctGeneID.mat'],'model')
exportToExcelFormat(model, [root '/scrap/SoyGEM_correctGeneID.xlsx']);

% Because of the large size of the model (particularly the large number of
% genes), it is probably more convenient to store intermediate model
% versions in YAML file format. This is more compact than SBML XML, and
% also has the benefit of being more easily readable (if you would want
% that).
writeYAMLmodel(model, [data '/templateModels/SoyGEM_correctGeneID.yml']);
