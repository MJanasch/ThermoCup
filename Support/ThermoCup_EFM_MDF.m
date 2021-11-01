%% Apply MDF to EFMs

function [MDF,model] = ThermoCup_EFM_MDF(model,metabolites,EFM)

%% Prune Model 
    Nr_Rxns = length(EFM);
    RunVar = Nr_Rxns:-1:1;
    for i = 1:Nr_Rxns
        if EFM(RunVar(i)) == 0 || isnan(model.reactions.dG0(RunVar(i)))
            
            model.rawS(:,RunVar(i))  = [];
            model.reactions.dG0(RunVar(i))  = [];
            model.reactions.rxns(RunVar(i))  = [];
            model.reactions.equations(RunVar(i))  = [];
            EFM(RunVar(i)) = [];
            
        end
    end

    EFM_Old = EFM;
    
    for j = 1:length(EFM)
        if EFM(j) < 0
            model.rawS(:,j) = model.rawS(:,j)*-1;
            EFM(j) = EFM(j)*-1;
            model.reactions.dG0(j) = model.reactions.dG0(j)*-1;
            model.reactions.equations(j) = strrep(model.reactions.equations(j),'<=>','<--');
        else
            model.reactions.equations(j) = strrep(model.reactions.equations(j),'<=>','-->');
        end
    end
    Rxns = model.reactions.rxns;
    model.metabolites = metabolites;

    
%% Calculate MDF for pathway

% Parameters
    T = 303.15;                     % Temperature [K]
    R = 0.008314;                   % Universal gas constant [kJ/(mol*K)]
    RT=R*T;

    [n,m] = size(model.rawS);

% Get standard delta G's (transpose to get from column vector to row vector
    dG0trans = transpose(model.reactions.dG0(:,1));

    solutionLP = MJanasch_ThermCup_SolveLPMDF(model,RT,n,m);%,beq);

    A = full(solutionLP.A);
    b = solutionLP.b;
    c = solutionLP.c;

    dG = model.reactions.dG0+RT*transpose(model.rawS)*solutionLP.z(1:n);
    dGtrans = transpose(dG);

    conc = exp(solutionLP.z(1:n));
    MDF = solutionLP.z(end,1)*RT;

    model.reactions.dG = dG;
    model.reactions.EFM = EFM;
    model.reactions.EFM_Old = EFM_Old;
%% Plot
    %disp(['Max-min Driving Force [kJ/mol] : ',num2str(round(MDF,4))])
end