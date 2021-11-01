%% Cupriavidus Core Model FBA
% Needs Cobra toolbox, check https://opencobra.github.io/cobratoolbox/stable/


%% Initialize CORBA toolbox
initCobraToolbox

%% Load core model
modelFile = '../Model/Cupriavidus_Core_ALL_FBA.xlsx';
CupriaCore_FBA = ThermoCup_ReadExcel_FBA(modelFile);

Nr_Rxns = length(CupriaCore_FBA.rxns);
%% Set constraints
% Fructose
%CupriaCore_FBA = changeRxnBounds(CupriaCore_FBA,'Ex_frc',-1,'b');

% Succinate
%CupriaCore_FBA = changeRxnBounds(CupriaCore_FBA,'Ex_suc',-1,'b');

% Formate
CupriaCore_FBA = changeRxnBounds(CupriaCore_FBA,'Ex_for',-1,'l');


%% Set to WT
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA,'XFPK1',0,'b');
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'XFPK2',0,'b');
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'XFPK3',0,'b');
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'aclAB',0,'b');

%% Set XFPK
CupriaCore_FBA_xfpk = changeRxnBounds(CupriaCore_FBA,'XFPK3',0,'b');
CupriaCore_FBA_xfpk = changeRxnBounds(CupriaCore_FBA_xfpk,'aclAB',0,'b');


%% Set aclAB
CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA,'XFPK1',0,'b');
CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA_aclAB,'XFPK2',0,'b');
CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA_aclAB,'XFPK3',0,'b');

%% Set Objective
CupriaCore_FBA_wt       = changeObjective(CupriaCore_FBA_wt,'Ex_Biomass',1);
CupriaCore_FBA_xfpk     = changeObjective(CupriaCore_FBA_xfpk,'Ex_Biomass',1);
CupriaCore_FBA_aclAB    = changeObjective(CupriaCore_FBA_aclAB,'Ex_Biomass',1);

%% Set Solver
changeCobraSolver('gurobi','LP');

%% Perform Optimization
FBAsolution_wt          = optimizeCbModel(CupriaCore_FBA_wt,'max');
FBAsolution_xfpk        = optimizeCbModel(CupriaCore_FBA_xfpk,'max');
FBAsolution_aclAB       = optimizeCbModel(CupriaCore_FBA_aclAB,'max');

Biomass_Max_wt          = FBAsolution_wt.obj;
Biomass_Max_xfpk        = FBAsolution_xfpk.obj;
Biomass_Max_aclAB       = FBAsolution_aclAB.obj;

CupriaCore_FBA_wt       = changeObjective(CupriaCore_FBA_wt,'Ex_PHB',1);
CupriaCore_FBA_xfpk     = changeObjective(CupriaCore_FBA_xfpk,'Ex_PHB',1);
CupriaCore_FBA_aclAB    = changeObjective(CupriaCore_FBA_aclAB,'Ex_PHB',1);

%% Get production envelope

% MolWeight = 0.18;     % Molecular weight of fructose
% MolWeight = 0.1181;   % Molecular weight of succinate
MolWeight = 0.045;      % Molecular weight of formate

x = linspace(100,0,1000);
for i = 1:length(x)
   CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'Ex_Biomass',Biomass_Max_wt*x(i)/100,'b');
   FBAsolution_wt = optimizeCbModel(CupriaCore_FBA_wt,'max');
   Result_wt(i,1) = Biomass_Max_wt*x(i)/(100*MolWeight);
   Result_wt(i,2) = FBAsolution_wt.obj/MolWeight;
   %Flux_Distribution(:,i) = FBAsolution.full;
   
   
   CupriaCore_FBA_xfpk = changeRxnBounds(CupriaCore_FBA_xfpk,'Ex_Biomass',Biomass_Max_xfpk*x(i)/100,'b');
   FBAsolution_xfpk = optimizeCbModel(CupriaCore_FBA_xfpk,'max');
   Result_xfpk(i,1) = Biomass_Max_xfpk*x(i)/(100*MolWeight);
   Result_xfpk(i,2) = FBAsolution_xfpk.obj/MolWeight;
   %Flux_Distribution(:,i) = FBAsolution.full;
   
   CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA_aclAB,'Ex_Biomass',Biomass_Max_aclAB*x(i)/100,'b');
   FBAsolution_aclAB = optimizeCbModel(CupriaCore_FBA_aclAB,'max');
   Result_aclAB(i,1) = Biomass_Max_aclAB*x(i)/(100*MolWeight);
   Result_aclAB(i,2) = FBAsolution_aclAB.obj/MolWeight;
   %Flux_Distribution(:,i) = FBAsolution.full;
end
