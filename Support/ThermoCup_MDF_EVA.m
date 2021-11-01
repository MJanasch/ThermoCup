%% Check EVA for concentrations

function Results_All = ThermoCup_MDF_EVA(modelOut,MDF,Leeway)
% loop through metabolites, maximize/minimize them
    model_EVA = modelOut;
    MDF_EVA_Orig = MDF;
    MDF_EVA_Leeway = MDF*Leeway; % Leeway = 0.9

    Nr_Met_EVA = length(model_EVA.metabolites.conc);

    for Met_Index = 1:Nr_Met_EVA
        % Maximize [Met]
        mima = 1; 
        
        % MDF at optimum
        [Results_Min_Orig] = FindExtreme_EVA(model_EVA,MDF_EVA_Orig,mima,Nr_Met_EVA,Met_Index);
        
        % MDF with margin
        [Results_Min_Leeway] = FindExtreme_EVA(model_EVA,MDF_EVA_Leeway,mima,Nr_Met_EVA,Met_Index);
        
        Results_Min.Conc_Out = min(Results_Min_Orig,Results_Min_Leeway);
        
        
        
        % Minimize [Met]
        mima = -1;
        
        % MDF at optimum
        [Results_Max_Orig] = FindExtreme_EVA(model_EVA,MDF_EVA_Orig,mima,Nr_Met_EVA,Met_Index);
        
        % MDF at optimum
        [Results_Max_Leeway] = FindExtreme_EVA(model_EVA,MDF_EVA_Leeway,mima,Nr_Met_EVA,Met_Index);
        
        
        Results_Max.Conc_Out = max(Results_Max_Orig,Results_Max_Leeway);
        
        
        Results_All(Met_Index,1) = Results_Min;
        Results_All(Met_Index,2) = Results_Max;
    end

    
    

end

function [Conc_Out] = FindExtreme_EVA(model_EVA,MDF_EVA,mima,Nr_Met_EVA,Met_Index)
    
    T = 303.15;                     % Temperature [K]
    R = 0.008314;                   % Universal gas constant [kJ/(mol*K)]
    RT=R*T;
    
    [n,m] = size(model_EVA.rawS);

%% Remove exchange reactions, as they are not part of the MDF algorithm
    S_no_Ex         = model_EVA.rawS;
    Strnsp_no_Ex    = model_EVA.Strnsp;
    rxns_dG         = model_EVA.reactions.rxns;

    
    RunVar = m:-1:1;
    for i = 1:length(RunVar)
        if strncmp(model_EVA.reactions.rxns(RunVar(i)),'Ex_',3)
            S_no_Ex(:,RunVar(i))            = [];
            Strnsp_no_Ex(:,RunVar(i))       = [];
            rxns_dG(RunVar(i))              = [];
        elseif strcmp(model_EVA.reactions.rxns(RunVar(i)),'BIOMASS')
            S_no_Ex(:,RunVar(i))            = [];
            Strnsp_no_Ex(:,RunVar(i))       = [];
            rxns_dG(RunVar(i))              = [];
        elseif strcmp(model_EVA.reactions.rxns(RunVar(i)),'Trp_H2O')
            S_no_Ex(:,RunVar(i))            = [];
            Strnsp_no_Ex(:,RunVar(i))       = [];
            rxns_dG(RunVar(i))              = [];
        elseif strcmp(model_EVA.reactions.rxns(RunVar(i)),'PhaCE')
            S_no_Ex(:,RunVar(i))            = [];
            Strnsp_no_Ex(:,RunVar(i))       = [];
            rxns_dG(RunVar(i))              = [];
        end
    end

    % Nr of transport reactions
    [Nr_Trxn_EVA,~] = size(Strnsp_no_Ex);
    
    % Nr of reactions without exchange reactions
    [~,m_no_Ex] = size(S_no_Ex);
    
%% Define conc vector to min/max
    Conc_Vector = zeros(Nr_Met_EVA,1);
    
    if mima == 1
        % Set value in c to 1
        Conc_Vector(Met_Index,1) = 1;

        c = [Conc_Vector;
            0;
            zeros(Nr_Met_EVA,1);
            zeros(Nr_Trxn_EVA,1)];
        
    elseif mima == -1
        % Set value in c to -1
        Conc_Vector(Met_Index,1) = -1;
        
        c = [Conc_Vector;
            0;
            zeros(Nr_Met_EVA,1);
            zeros(Nr_Trxn_EVA,1)];
    else
        disp('minmax needs to either 1 or -1 to find the maximum or minimum of the concentration range, respectively.')
    end
    
    %% Define A and b according to ThermoCup_MDF
    A = model_EVA.optimized.A;
%     A = [transpose(S_no_Ex)         ones(m_no_Ex,1)             transpose(S_no_Ex)          transpose(Strnsp_no_Ex);
%          eye(n)                     zeros(n,1)                  zeros(n,n)                  zeros(n,Nr_Trxn_EVA_test);
%         -eye(n)                     zeros(n,1)                  zeros(n,n)                  zeros(n,Nr_Trxn_EVA_test);
%         zeros(n,n)                  zeros(n,1)                  eye(n)                      zeros(n,Nr_Trxn_EVA_test);
%         zeros(n,n)                  zeros(n,1)                  -eye(n)                     zeros(n,Nr_Trxn_EVA_test);
%         zeros(Nr_Trxn_EVA_test,n)   zeros(Nr_Trxn_EVA_test,1)	zeros(Nr_Trxn_EVA_test,n)  	eye(Nr_Trxn_EVA_test);
%         zeros(Nr_Trxn_EVA_test,n)   zeros(Nr_Trxn_EVA_test,1)	zeros(Nr_Trxn_EVA_test,n)	-eye(Nr_Trxn_EVA_test)];

    b = model_EVA.optimized.b;
%     b = [zeros(m_no_Ex,1);
%          log(uconc);
%         -log(lconc);
%          dfG/RT;
%         -dfG/RT;
%          dtG/RT;
%         -dtG/RT];
%     
    
    
    %% Set Aeq and beq so that B = MDF
    Aeq = [zeros(1,Nr_Met_EVA) 1 zeros(1,Nr_Met_EVA) zeros(1,Nr_Trxn_EVA)];
   
    % Define beq
    beq = MDF_EVA/RT;
    
    % Solve
    options = optimoptions('linprog','Display','none');
    X = linprog(transpose(c),A,b,Aeq,beq,[],[],options);
    Conc_Out = 1000*exp(X(Met_Index)); % Multiply by 1000 to get from M to mM
    
end