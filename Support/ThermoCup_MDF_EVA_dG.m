%% Check EVA for dG

function Results_All = ThermoCup_MDF_EVA_dG(modelOut,MDF,Leeway)
    model_EVA = modelOut;
    MDF_EVA_Orig = MDF;
    MDF_EVA_Leeway = MDF*Leeway; % Leeway = 0.9 reducing MDF by 10% to account for uncertainty and potentially releave zero-range metabolite concentration constraints

    
    % Remove exchange reactions, as they are not part of the MDF algorithm
    
    [n,m] = size(model_EVA.rawS);
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


    [~,m_no_Ex] = size(S_no_Ex);
    
    % Loop through reactions, maximize/minimize them
    for Rxn_Index = 1:m_no_Ex
        % Maximize [dG]
        mima = 1; 
    
        % MDF at optimum
        [Results_Min_Orig] = FindExtreme_EVA_dG(model_EVA,MDF_EVA_Orig,mima,S_no_Ex,Strnsp_no_Ex,Rxn_Index);
        
        % MDF with margin
        [Results_Min_Leeway] = FindExtreme_EVA_dG(model_EVA,MDF_EVA_Leeway,mima,S_no_Ex,Strnsp_no_Ex,Rxn_Index);
        
        
        Results_Min = min(Results_Min_Orig.dG,Results_Min_Leeway.dG);
        
        
        % Minimize [dG]
        mima = -1;
        
        % MDF at optimum
        [Results_Max_Orig] = FindExtreme_EVA_dG(model_EVA,MDF_EVA_Orig,mima,S_no_Ex,Strnsp_no_Ex,Rxn_Index);
        
        % MDF with margin
        [Results_Max_Leeway] = FindExtreme_EVA_dG(model_EVA,MDF_EVA_Leeway,mima,S_no_Ex,Strnsp_no_Ex,Rxn_Index);
        
        Results_Max = max(Results_Max_Orig.dG,Results_Max_Leeway.dG);
        
        
        
        Results_All(Rxn_Index,1) = Results_Min;
        Results_All(Rxn_Index,2) = Results_Max;
    end
end

function [Results_Range] = FindExtreme_EVA_dG(model_EVA,MDF_EVA,mima,S_no_Ex,Strnsp_no_Ex,Rxn_Index)
    
    T = 303.15;                     % Temperature [K]
    R = 0.008314;                   % Universal gas constant [kJ/(mol*K)]
    RT=R*T;
    
    [n,m] = size(model_EVA.rawS);

%% Get physiological boundaries
    dfG     = model_EVA.metabolites.dfG;

%     lconc   = model_EVA.metabolites.lconc;
%     uconc   = model_EVA.metabolites.uconc;

    dtG     = model_EVA.reactions.dtG;

    % Nr of transport reactions
    [Nr_Trxn_EVA_test,~] = size(Strnsp_no_Ex);
    
    % Nr of reactions without exchange reactions
    [n_no_Ex,m_no_Ex] = size(S_no_Ex);
    
%% Define vector to min/max
    Obj_Vector = S_no_Ex(:,Rxn_Index);
    Obj_Vector_Trnsp = Strnsp_no_Ex(:,Rxn_Index);
    
    c = [Obj_Vector;
        0;
        Obj_Vector;
        Obj_Vector_Trnsp];

    %% Define A and b according to MDF 2.0
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

    Aeq = [zeros(1,n_no_Ex) 1 zeros(1,n_no_Ex) zeros(1,Nr_Trxn_EVA_test)];

    % Define beq
    beq = MDF_EVA/RT;

    % Solve
    options = optimoptions('linprog','Display','none');
    if mima == 1
        X = linprog(transpose(-c),A,b,Aeq,beq,[],[],options);
        
    elseif mima == -1
        X = linprog(transpose(c),A,b,Aeq,beq,[],[],options);
    else
        disp('minmax needs to either 1 or -1 to find the maximum or minimum of the concentration range, respectively.')
    end
    
    Results_Range.dG = transpose(S_no_Ex(:,Rxn_Index))*dfG+RT*transpose(S_no_Ex(:,Rxn_Index))*X(1:n)+transpose(Strnsp_no_Ex(:,Rxn_Index))*dtG;
   
end