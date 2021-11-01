%% MDF for small networks

%% Load model and data

addpath('Model/')
modelFile = '/Model/Cupriavidus_Core_XII_PHB.xlsx';
CupriaCore = ThermoCup_ReadExcel(modelFile);

addpath('Data/')
load('ThermoCup_All_Enumerated_EFMs.mat');

addpath('Support/')


Nr_Rxns = length(mnet_For.reactionNames);

%% Create EFM for canonical TCA
Rxn_Names_CanonTCA = {'PYK';
                    'PDH';
                    'PYRC';
                    'SUCD';
                    'FUMC';
                    'MDH';
                    'CITS';
                    'ACH'};

EFM_CanonTCA = zeros(Nr_Rxns,1);

EFM_CanonTCA(ThermoCup_FindIndex(mnet_For.reactionNames,Rxn_Names_CanonTCA)) = 1;

%% Create EFM for MDH reversal
Rxn_Names_revMDH_fwd = {'PYK';
                    'PDH';
                    'SUCD';
                    'FUMC';
                    'CITS';
                    'ACH';
                    'MAE'};

Rxn_Names_revMDH_rev = {'PYRC';
                    'MDH'};
           
EFM_revMDH = zeros(Nr_Rxns,1);
EFM_revMDH(ThermoCup_FindIndex(mnet_For.reactionNames,Rxn_Names_revMDH_fwd)) = 1;
EFM_revMDH(ThermoCup_FindIndex(mnet_For.reactionNames,Rxn_Names_revMDH_rev)) = -1;

%% Calculate MDF and dGs
[MDF_canonTCA,modelOut_canonTCA] = ThermoCup_MDF_SmallNetwork(CupriaCore.Variant_For,CupriaCore.metabolites,CupriaCore.cofaratios,EFM_CanonTCA);
[MDF_revMDH,modelOut_revMDH] = ThermoCup_MDF_SmallNetwork(CupriaCore.Variant_For,CupriaCore.metabolites,CupriaCore.cofaratios,EFM_revMDH);

EVA_dG.Results_canonTCA = ThermoCup_MDF_EVA_dG(modelOut_canonTCA,MDF_canonTCA,0.999);
EVA_dG.Results_revMDH = ThermoCup_MDF_EVA_dG(modelOut_revMDH,MDF_revMDH,0.999);

%% dG ranges in better form
for i = 1:length(EVA_dG.Results_canonTCA)
    DF_canonTCA(i,1) = -round(EVA_dG.Results_canonTCA(i,1),4);
    DF_canonTCA(i,2) = -round(EVA_dG.Results_canonTCA(i,2),4);
end

for i = 1:length(EVA_dG.Results_revMDH)
    DF_revMDH(i,1) = -round(EVA_dG.Results_revMDH(i,1),4);
    DF_revMDH(i,2) = -round(EVA_dG.Results_revMDH(i,2),4);
end