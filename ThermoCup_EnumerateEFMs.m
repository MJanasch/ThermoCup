%% Main file


%% Load model
addpath('Model/')
addpath('Support/')
modelFile = '/Model/Cupriavidus_Core_XII_PHB.xlsx';
CupriaCore = ThermoCup_ReadExcel(modelFile);

%% Prepare structure for EFM enumeration via EFMtool

% Need to move to EFMtool directory or add path
% cd('.../efmtool/')
% stru.stoich             = CupriaCore.Variant_For.rawS;                  % Change to CupriaCore.Variant_Frc.rawS or CupriaCore.Variant_Suc.rawS respecetively
% stru.reversibilities    = CupriaCore.Variant_For.reactions.reversible;  % Change to CupriaCore.Variant_Frc.reactions.reversible or CupriaCore.Variant_Suc.reactions.reversible respecetively
% stru.metaboliteNames    = CupriaCore.metabolites.mets;                  % Can stay the same for all 3 variants
% stru.reactionNames      = CupriaCore.Variant_For.reactions.rxns;        % Change to CupriaCore.Variant_Frc.reactions.rxns or CupriaCore.Variant_Suc.reactions.rxns, respectively
% 
% % Enumerate EFMs via EFMtool
% tic
% mnet_For = CalculateFluxModes(stru);                                    % Change to mnet_Frc or mnet_Suc, repsectively
% toc

%% Load enumerated EFMs
% everything under "Prepare structure for EFM enumeration via EFMtool" has been performed and all three mnet-variants have been saved
addpath('Data/')
load('ThermoCup_All_Enumerated_EFMs.mat');

%% Check for biomass and PHB production
% Choose the substrate
% Substrate = 'Frc';
% Substrate = 'Suc';
Substrate = 'For';


if exist('Substrate','var') == 1
    if strcmp(Substrate,'Frc')
        mnet = mnet_Frc;
        MolWeight = 0.18;     % Molecular weight of fructose
        Xfpk1 = 72;
        Xfpk2 = 73;
        AclAB = 74;
        Model_Variant = CupriaCore.Variant_Frc;
        disp('Substrate is Frc')
    elseif strcmp(Substrate,'Suc')
        mnet = mnet_Suc;
        MolWeight = 0.1181;   % Molecular weight of succinate
        Xfpk1 = 72;
        Xfpk2 = 73;
        AclAB = 74;
        Model_Variant = CupriaCore.Variant_Suc;
        disp('Substrate is Suc')
    elseif strcmp(Substrate,'For')
        mnet = mnet_For;
        MolWeight = 0.045;      % Molecular weight of formate
        Xfpk1 = 73;
        Xfpk2 = 74;
        AclAB = 75;
        Model_Variant = CupriaCore.Variant_For;
        disp('Substrate is For')
    else
        disp('Choose a right substrate')
    end
else
    disp('Define a substrate, please choose between "Frc", "Suc" and "For".')
end

Nr_Rxns = length(mnet.reactionNames);


%% Select EFMs for reaction usage
% mnet.efms(4,i) is always biomass
% mnet.efms(9,i) is always PHB exchange

% for fructose and succinate
    % mnet.efms(Xfpk1,i) is XFPK1
    % mnet.efms(Xfpk2,i) is XFPK2
    % mnet.efms(AclAB,i) is aclAB

% for formate
    % mnet.efms(Xfpk2,i) is XFPK1
    % mnet.efms(AclAB,i) is XFPK2
    % mnet.efms(75,i) is aclAB
 
%% Wildtype H16
tic
r=0;
r_yes = 0;
for i = 1:length(mnet.efms)
    %if mnet.efms(9,i) ~= 0 && mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) == 0 && mnet.efms(AclAB,i) == 0
    if (((mnet.efms(4,i) ~= 0) && mnet.efms(9,i) ~= 0) || (mnet.efms(4,i) ~= 0 && mnet.efms(9,i) == 0) || (mnet.efms(4,i) == 0 && mnet.efms(9,i) ~= 0)) && mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) == 0 && mnet.efms(AclAB,i) == 0
        r=r+1;
        Result(i) = 1;
        Index(r) = i;
        [MDF(r)] = ThermoCup_MDF(Model_Variant,CupriaCore.metabolites,CupriaCore.cofaratios,mnet.efms(:,i));
        
        if MDF(r) > 0    % If the EFM is thermodynamically feasible
            r_yes = r_yes +1;
            MDF_yes(r_yes) = MDF(r);
            % Perform Energy Variability Analysis (EVA)
            %EVA(r_yes).Results = ThermoCup_MDF_EVA(modelOut(r),MDF(r),0.9); % Check ranges for metabolite concentrations
            %EVA_dG(r_yes).Results = ThermoCup_MDF_EVA_dG(modelOut(r),MDF(r),0.9); % Check ranges for dG's
             
            for o = 1:Nr_Rxns
                if mnet.efms(o,i) < 0
                    EFM_Count(o,r_yes) = -1;
                elseif mnet.efms(o,i) > 0
                    EFM_Count(o,r_yes) = 1;
                else
                    EFM_Count(o,r_yes) = 0;
                end
            end
            
            yield_Biomass(r_yes)    = mnet.efms(4,i)/(mnet.efms(1,i)*MolWeight); 
            yield_PHB(r_yes)        = mnet.efms(9,i)/(mnet.efms(1,i)*MolWeight);
            N_Uptake(r_yes)         = mnet.efms(8,i);
            if rem(r_yes,1000) == 0
                disp(['Number of feasible EFMs found: ',num2str(r_yes)])
                toc
            end
        end
    else
        Result(i) = 0;
    end
    
end
disp(['Total number of feasible EFMs found: ',num2str(r_yes)])
toc

%% H16-AclAB
w=0;
w_yes = 0;
for i = 1:length(mnet.efms)
    %if mnet.efms(9,i) ~= 0 && mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) == 0 && mnet.efms(AclAB,i) ~= 0
    if (((mnet.efms(4,i) ~= 0) && mnet.efms(9,i) ~= 0) || (mnet.efms(4,i) ~= 0 && mnet.efms(9,i) == 0) || (mnet.efms(4,i) == 0 && mnet.efms(9,i) ~= 0)) && mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) == 0 && mnet.efms(AclAB,i) ~= 0
        w=w+1;
        Result_aclAB(i) = 1;
        Index_aclAB(w) = i;
        [MDF_aclAB(w)] = ThermoCup_MDF(Model_Variant,CupriaCore.metabolites,CupriaCore.cofaratios,mnet.efms(:,i));
        if MDF_aclAB(w) > 0     % If the EFM is thermodynamically feasible
            w_yes = w_yes + 1;
            MDF_aclAB_yes(w_yes) = MDF_aclAB(w);
            %EVA_aclAB(w).Results = ThermoCup_MDF_EVA(modelOut_aclAB(w),MDF_aclAB(w),0.9);
            %EVA_dG_aclAB(w).Results = ThermoCup_MDF_EVA_dG(modelOut_aclAB(w),MDF_aclAB(w),0.9);
            
            for o = 1:Nr_Rxns
                if mnet.efms(o,i) < 0
                    EFM_Count_aclAB(o,w_yes) = -1;
                elseif mnet.efms(o,i) > 0
                    EFM_Count_aclAB(o,w_yes) = 1;
                else
                    EFM_Count_aclAB(o,w_yes) = 0;
                end
            end
            yield_Biomass_aclAB(w_yes)      = mnet.efms(4,i)/(mnet.efms(1,i)*MolWeight);
            yield_PHB_aclAB(w_yes)          = mnet.efms(9,i)/(mnet.efms(1,i)*MolWeight);
            N_Uptake_aclAB(w_yes)           = mnet.efms(8,i);
        
            if rem(w_yes,1000) == 0
                disp(['Number of feasible AclAB-EFMs found: ',num2str(w_yes)])
                toc
            end
        
        end
        
        

        
    else
        Result_aclAB(i) = 0;
    end 
end
toc

%% H16-Xfpk
v=0;
v_yes=0;
for i = 1:length(mnet.efms)
    %if mnet.efms(9,i) ~= 0 && (((mnet.efms(Xfpk1,i) ~= 0) && mnet.efms(Xfpk2,i) ~= 0) || (mnet.efms(Xfpk1,i) ~= 0 && mnet.efms(Xfpk2,i) == 0) || (mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) ~= 0)) && mnet.efms(AclAB,i) == 0
    if (((mnet.efms(4,i) ~= 0) && mnet.efms(9,i) ~= 0) || (mnet.efms(4,i) ~= 0 && mnet.efms(9,i) == 0) || (mnet.efms(4,i) == 0 && mnet.efms(9,i) ~= 0)) && (((mnet.efms(Xfpk1,i) ~= 0) && mnet.efms(Xfpk2,i) ~= 0) || (mnet.efms(Xfpk1,i) ~= 0 && mnet.efms(Xfpk2,i) == 0) || (mnet.efms(Xfpk1,i) == 0 && mnet.efms(Xfpk2,i) ~= 0)) && mnet.efms(AclAB,i) == 0
        v=v+1;
        Result_xfpk(i) = 1;
        Index_xfpk(v) = i;
       [MDF_xfpk(v)] = ThermoCup_MDF(Model_Variant,CupriaCore.metabolites,CupriaCore.cofaratios,mnet.efms(:,i));
        if MDF_xfpk(v) > 0     % If the EFM is thermodynamically feasible
            v_yes = v_yes +1;
            MDF_xfpk_yes(v_yes) = MDF_xfpk(v);
            %EVA_xfpk(v_yes).Results = ThermoCup_MDF_EVA(modelOut_xfpk(v),MDF_xfpk(v),0.9);
            %EVA_dG_xfpk(v_yes).Results = ThermoCup_MDF_EVA_dG(modelOut_xfpk(v),MDF_xfpk(v),0.9);
            for o = 1:Nr_Rxns
                if mnet.efms(o,i) < 0
                    EFM_Count_xfpk(o,v_yes) = -1;
                elseif mnet.efms(o,i) > 0
                    EFM_Count_xfpk(o,v_yes) = 1;
                else
                    EFM_Count_xfpk(o,v_yes) = 0;
                end
            end
            yield_Biomass_xfpk(v_yes)      = mnet.efms(4,i)/(mnet.efms(1,i)*MolWeight);
            yield_PHB_xfpk(v_yes)          = mnet.efms(9,i)/(mnet.efms(1,i)*MolWeight);
            N_Uptake_xfpk(v_yes)           = mnet.efms(8,i);
            if rem(v_yes,1000) == 0
                disp(['Number of feasible Xfpk-EFMs found: ',num2str(v_yes)])
                toc
            end
        end 
        
        

        
    else
        Result_xfpk(i) = 0;
    end 
end
toc
