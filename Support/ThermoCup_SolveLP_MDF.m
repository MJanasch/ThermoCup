function [solution,S_no_Ex,Strnsp_no_Ex,rxns_dG] = ThermoCup_SolveLP_MDF(modelPruned,RT,dtG)

[n,m] = size(modelPruned.rawS);
% n: Nr of metabolites
% m: Nr of reactions before removing
%% Remove exchange reactions, as they are not part of the MDF algorithm, 
% and remove biomass equation (assumed to be sufficiently forward driven under all conditions
% and remove water transport
% and remove PHB synthase as there is no effective concentration of PHB polymers

S_no_Ex         = modelPruned.rawS;
Strnsp_no_Ex    = modelPruned.Strnsp;
rxns_dG         = modelPruned.reactions.rxns;

RunVar = m:-1:1;
for i = 1:length(RunVar)
    if strncmp(modelPruned.reactions.rxns(RunVar(i)),'Ex_',3)
        S_no_Ex(:,RunVar(i))            = [];
        Strnsp_no_Ex(:,RunVar(i))       = [];
        rxns_dG(RunVar(i))              = [];
    elseif strcmp(modelPruned.reactions.rxns(RunVar(i)),'BIOMASS')
        S_no_Ex(:,RunVar(i))            = [];
        Strnsp_no_Ex(:,RunVar(i))       = [];
        rxns_dG(RunVar(i))              = [];
    elseif strncmp(modelPruned.reactions.rxns(RunVar(i)),'Trp_H2O',7)
        S_no_Ex(:,RunVar(i))            = [];
        Strnsp_no_Ex(:,RunVar(i))       = [];
        rxns_dG(RunVar(i))              = [];
    elseif strcmp(modelPruned.reactions.rxns(RunVar(i)),'PhaCE')
        S_no_Ex(:,RunVar(i))            = [];
        Strnsp_no_Ex(:,RunVar(i))       = [];
        rxns_dG(RunVar(i))              = [];
    end
end

[~,m_no_Ex] = size(S_no_Ex);
% m_no_Ex: Nr of reactions minus the ones removed above

%% Number of transport reactions
Nr_Trnsp = length(dtG);

%% Create S-matrix for ratio constraints
% S_ratios has to be Nr_Rations x Nr_metabolites, here 4xn
Nr_Ratios = length(modelPruned.cofaratios.ratios);

S_ratios = zeros(n,Nr_Ratios);
for j = 1:Nr_Ratios
    Ratio_Split = split(modelPruned.cofaratios.ratios(j),"/");
    Over = Ratio_Split(1);
    Under = Ratio_Split(2);
    S_ratios(MJanasch_FindIndex(modelPruned.metabolites.mets,Over),j) = 1;
    S_ratios(MJanasch_FindIndex(modelPruned.metabolites.mets,Under),j) = -1;
end
    
%% Get thermodynamic data
dfG     = modelPruned.metabolites.dfG;

lconc   = modelPruned.metabolites.lconc;
uconc   = modelPruned.metabolites.uconc;

lratio = modelPruned.cofaratios.lratio;
uratio = modelPruned.cofaratios.uratio;

%% Formulate LP Problem
% Define Objective
c = [zeros(n,1);
     1;
     zeros(n,1);
     zeros(Nr_Trnsp,1)];

% Define A
A = [transpose(S_no_Ex)     ones(m_no_Ex,1)     transpose(S_no_Ex)  transpose(Strnsp_no_Ex);
     eye(n)                 zeros(n,1)          zeros(n,n)          zeros(n,Nr_Trnsp);
    -eye(n)                 zeros(n,1)          zeros(n,n)          zeros(n,Nr_Trnsp);
    zeros(n,n)              zeros(n,1)          eye(n)              zeros(n,Nr_Trnsp);
    zeros(n,n)              zeros(n,1)          -eye(n)             zeros(n,Nr_Trnsp);
    zeros(Nr_Trnsp,n)       zeros(Nr_Trnsp,1)	zeros(Nr_Trnsp,n)  	eye(Nr_Trnsp);
    zeros(Nr_Trnsp,n)       zeros(Nr_Trnsp,1)	zeros(Nr_Trnsp,n)	-eye(Nr_Trnsp);
    transpose(S_ratios)     zeros(Nr_Ratios,1)  zeros(Nr_Ratios,n)  zeros(Nr_Ratios,Nr_Trnsp);
    transpose(-S_ratios)    zeros(Nr_Ratios,1)  zeros(Nr_Ratios,n)  zeros(Nr_Ratios,Nr_Trnsp)];

% Define b
b = [zeros(m_no_Ex,1);
     log(uconc);
    -log(lconc);
    dfG/RT;
    -dfG/RT;
    dtG/RT;
    -dtG/RT;
    log(uratio);
    -log(lratio)];

solution.c=c;
solution.A=A;
solution.b=b;

options = optimoptions('linprog','Display','none');
% Primal Problem: Finding MDF
X = linprog(transpose(-c),A,b,[],[],[],[],options);
solution.z=X;

end