%% Apply ERA-MDF to EFMs

function [MDF,modelPruned] = ThermoCup_MDF(model,metabolites,cofaratios,EFM)

%% Parameters
    T = 303.15;                     % Temperature [K]
    R = 0.008314;                   % Universal gas constant [kJ/(mol*K)]
    RT=R*T;

    pH_in   = 7.5;
    pH_out  = 7;
    
    membrane_potential = -0.15; % [V]

    F = 96.48533212; % [C/mmol]

%% Prune Model, i.e. remove all zero flux reactions and all unused metabolites
    Nr_Rxns = length(EFM);
    modelPruned = model;
    EFM_Old = EFM;
    EFM_Pruned = EFM;
    RunVar = Nr_Rxns:-1:1;
    
    
    for i = 1:Nr_Rxns
        if EFM_Pruned(RunVar(i)) == 0
            modelPruned.rawS(:,RunVar(i))               = [];
            modelPruned.reactions.rxns(RunVar(i))       = [];
            modelPruned.reactions.equations(RunVar(i))  = [];
            modelPruned.reactions.transport(RunVar(i))  = [];
            EFM_Pruned(RunVar(i))                       = [];
        end
    end
    
    modelPruned.metabolites = metabolites;
    
    Nr_Metabolites = length(modelPruned.metabolites.mets);
    RunVar_Metabolites = Nr_Metabolites:-1:1;
    for j = 1:Nr_Metabolites
        if all(modelPruned.rawS(RunVar_Metabolites(j),:) == 0)
            modelPruned.rawS(RunVar_Metabolites(j),:)                       = [];
            modelPruned.metabolites.metNames(RunVar_Metabolites(j))         = [];
            modelPruned.metabolites.mets(RunVar_Metabolites(j))             = [];
            modelPruned.metabolites.dfG(RunVar_Metabolites(j))              = [];
            modelPruned.metabolites.dfGUncertainty(RunVar_Metabolites(j))   = [];
            modelPruned.metabolites.lconc(RunVar_Metabolites(j))            = [];
            modelPruned.metabolites.uconc(RunVar_Metabolites(j))            = [];
        end
    end
        
%% Add cofactor ratios to model
    modelPruned.cofaratios = cofaratios;
    
%% Add transport energies to transport reactions

[dtG_Names,dtG_Values] = ThermoCup_CalculateTransportEnergy(RT,pH_in,pH_out,membrane_potential,F);

[n,m] = size(modelPruned.rawS);

Nr_trnsp = sum(modelPruned.reactions.transport);

modelPruned.Strnsp = zeros(Nr_trnsp,m);

Trnsp_RunVar=0;
for i = 1:length(modelPruned.reactions.rxns)
    for j = 1:length(dtG_Names)
        if strcmp(modelPruned.reactions.rxns{i},dtG_Names{j})
            Trnsp_RunVar = Trnsp_RunVar+1;
            modelPruned.Strnsp(Trnsp_RunVar,i) = 1;
            dtG(Trnsp_RunVar,1) = dtG_Values(j);
        end
    end
end


%% Make all reactions forward (no negative fluxes)    

    for j = 1:length(EFM_Pruned)
        if EFM_Pruned(j) < 0
            modelPruned.rawS(:,j) = modelPruned.rawS(:,j)*-1;
            EFM_Pruned(j) = EFM_Pruned(j)*-1;
            modelPruned.Strnsp(:,j) = modelPruned.Strnsp(:,j)*-1;
            modelPruned.reactions.equations(j) = strrep(modelPruned.reactions.equations(j),'<=>','<--');
            modelPruned.reactions.rxns{j} = [modelPruned.reactions.rxns{j},'_rev'];
        else
            modelPruned.reactions.equations(j) = strrep(modelPruned.reactions.equations(j),'<=>','-->');
        end
    end

%% Calculate MDF for pathway

	[n,m] = size(modelPruned.rawS);

    [solutionLP,S_MOD,Strnsp_MOD,rxns_dG] = ThermoCup_SolveLP_MDF(modelPruned,RT,dtG);
    
    conc = exp(solutionLP.z(1:n));
    MDF = solutionLP.z(n+1,1)*RT;

    modelPruned.metabolites.conc        = conc;
    modelPruned.reactions.dtG           = dtG;
    modelPruned.S_MOD                   = S_MOD;
    modelPruned.Strnsp_MOD              = Strnsp_MOD;
    modelPruned.reactions.rxns_dG       = rxns_dG;
    modelPruned.optimized               = solutionLP;

end