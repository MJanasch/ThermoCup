%% Cupriavidus Core Model FBA
% Needs Cobra toolbox, check https://opencobra.github.io/cobratoolbox/stable/


%% Initialize CORBA toolbox
% initCobraToolbox

%% Load core model
modelFile = '../Model/Cupriavidus_Core_ALL_FBA_II.xlsx';
CupriaCore_FBA = ThermoCup_ReadExcel_FBA(modelFile);

Nr_Rxns = length(CupriaCore_FBA.rxns);
%% Set constraints
% Fructose
CupriaCore_FBA = changeRxnBounds(CupriaCore_FBA,'Ex_frc',-1,'b');
MolWeight = 0.18;     % Molecular weight of fructose



% Formate
%CupriaCore_FBA = changeRxnBounds(CupriaCore_FBA,'Ex_for',-1,'l');

%MolWeight = 0.045;      % Molecular weight of formate

%% Set to WT
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA,'XFPK1',0,'b');
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'XFPK2',0,'b');
CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'aclAB',0,'b');

%% Set XFPK
CupriaCore_FBA_xfpk = changeRxnBounds(CupriaCore_FBA,'aclAB',0,'b');


%% Set aclAB
CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA,'XFPK1',0,'b');
CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA_aclAB,'XFPK2',0,'b');

%% Set Objective
CupriaCore_FBA_wt       = changeObjective(CupriaCore_FBA_wt,'Ex_Biomass',1);
CupriaCore_FBA_xfpk     = changeObjective(CupriaCore_FBA_xfpk,'Ex_Biomass',1);
CupriaCore_FBA_aclAB    = changeObjective(CupriaCore_FBA_aclAB,'Ex_Biomass',1);

%% Set Solver
changeCobraSolver('ibm_cplex','LP');

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

MolWeight_PHB = 0.086; % Molecular weight of PHB monomer

x = linspace(1,0,1000); % in 1000 steps, go down from 1 to 0
for i = 1:length(x)
   CupriaCore_FBA_wt = changeRxnBounds(CupriaCore_FBA_wt,'Ex_Biomass',Biomass_Max_wt*x(i),'b');
   FBAsolution_wt = optimizeCbModel(CupriaCore_FBA_wt,'max');
   Result_wt(i,1) = Biomass_Max_wt*x(i)/MolWeight;
   Result_wt(i,2) = FBAsolution_wt.obj*MolWeight_PHB/MolWeight;
   if Result_wt(i,2) < 0
       Result_wt(i,2) = 0;
   end
   %Flux_Distribution(:,i) = FBAsolution.full;
   
   
   CupriaCore_FBA_xfpk = changeRxnBounds(CupriaCore_FBA_xfpk,'Ex_Biomass',Biomass_Max_xfpk*x(i),'b');
   FBAsolution_xfpk = optimizeCbModel(CupriaCore_FBA_xfpk,'max');
   Result_xfpk(i,1) = Biomass_Max_xfpk*x(i)/MolWeight;
   Result_xfpk(i,2) = FBAsolution_xfpk.obj*MolWeight_PHB/MolWeight;
   if Result_xfpk(i,2) < 0
       Result_xfpk(i,2) = 0;
   end
   %Flux_Distribution(:,i) = FBAsolution.full;
   
   CupriaCore_FBA_aclAB = changeRxnBounds(CupriaCore_FBA_aclAB,'Ex_Biomass',Biomass_Max_aclAB*x(i),'b');
   FBAsolution_aclAB = optimizeCbModel(CupriaCore_FBA_aclAB,'max');
   Result_aclAB(i,1) = Biomass_Max_aclAB*x(i)/MolWeight;
   Result_aclAB(i,2) = FBAsolution_aclAB.obj*MolWeight_PHB/MolWeight;
   if Result_aclAB(i,2) < 0
       Result_aclAB(i,2) = 0;
   end
   
   %Flux_Distribution(:,i) = FBAsolution.full;
end


figure(1)
plot(Result_wt(:,1),Result_wt(:,2))
hold on
plot(Result_xfpk(:,1),Result_xfpk(:,2))
plot(Result_aclAB(:,1),Result_aclAB(:,2))

load handel
sound(y,Fs)