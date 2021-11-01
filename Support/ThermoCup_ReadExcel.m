
function model = ThermoCup_ReadExcel(fileName)
data_raw = importdata(fileName);

% Reactions
    model.Variant_Frc.reactions.rxns        = data_raw.textdata.Reactions_Frc(2:end,1);
    model.Variant_Frc.reactions.rxnNames    = data_raw.textdata.Reactions_Frc(2:end,3);
    model.Variant_Frc.reactions.equations   = data_raw.textdata.Reactions_Frc(2:end,4);
    model.Variant_Frc.reactions.lb          = [];
    model.Variant_Frc.reactions.ub          = [];
    model.Variant_Frc.reactions.eccodes     = data_raw.textdata.Reactions_Frc(2:end,6);
    model.Variant_Frc.reactions.dG0         = data_raw.data.Reactions_Frc(:,1);
    model.Variant_Frc.reactions.reversible  = data_raw.data.Reactions_Frc(:,3);
    model.Variant_Frc.reactions.transport   = data_raw.data.Reactions_Frc(:,4);
    
    model.Variant_Suc.reactions.rxns        = data_raw.textdata.Reactions_Suc(2:end,1);
    model.Variant_Suc.reactions.rxnNames    = data_raw.textdata.Reactions_Suc(2:end,3);
    model.Variant_Suc.reactions.equations   = data_raw.textdata.Reactions_Suc(2:end,4);
    model.Variant_Suc.reactions.lb          = [];
    model.Variant_Suc.reactions.ub          = [];
    model.Variant_Suc.reactions.eccodes     = data_raw.textdata.Reactions_Suc(2:end,6);
    model.Variant_Suc.reactions.dG0         = data_raw.data.Reactions_Suc(:,1);
    model.Variant_Suc.reactions.reversible  = data_raw.data.Reactions_Suc(:,3);
    model.Variant_Suc.reactions.transport   = data_raw.data.Reactions_Suc(:,4);
    
    model.Variant_For.reactions.rxns        = data_raw.textdata.Reactions_For(2:end,1);
    model.Variant_For.reactions.rxnNames    = data_raw.textdata.Reactions_For(2:end,3);
    model.Variant_For.reactions.equations   = data_raw.textdata.Reactions_For(2:end,4);
    model.Variant_For.reactions.lb          = [];
    model.Variant_For.reactions.ub          = [];
    model.Variant_For.reactions.eccodes     = data_raw.textdata.Reactions_For(2:end,6);
    model.Variant_For.reactions.dG0         = data_raw.data.Reactions_For(:,1);
    model.Variant_For.reactions.reversible  = data_raw.data.Reactions_For(:,3);
    model.Variant_For.reactions.transport   = data_raw.data.Reactions_For(:,4);
    
%% Metabolites

    Nr_R = 0;
    Nr_Mets = length(data_raw.textdata.Metabolites(2:end,1));
    metnames_raw = data_raw.textdata.Metabolites(2:end,1);
    %model.metabolites.mets=data_raw.textdata.Metabolites(2:end,1);
    for i = 1:Nr_Mets
        if ~contains(metnames_raw{i},'/')
            model.metabolites.mets(i,1)                  = data_raw.textdata.Metabolites(i+1,1);
            model.metabolites.metNames(i,1)              = data_raw.textdata.Metabolites(i+1,2);
            model.metabolites.lconc(i,1)                 = data_raw.data.Metabolites(i,1);%*model.annotation.defaultLCB;
            model.metabolites.uconc(i,1)                 = data_raw.data.Metabolites(i,2);%*model.annotation.defaultUCB;
            model.metabolites.dfG(i,1)                   = data_raw.data.Metabolites(i,3);
            model.metabolites.dfGUncertainty(i,1)        = data_raw.data.Metabolites(i,4);
        elseif contains(metnames_raw{i},'/')
            Nr_R = Nr_R + 1;
            model.cofaratios.ratios(Nr_R,1)                = data_raw.textdata.Metabolites(i+1,1);
            model.cofaratios.ratioNames(Nr_R,1)            = data_raw.textdata.Metabolites(i+1,2);
            model.cofaratios.lratio(Nr_R,1)                = data_raw.data.Metabolites(i,1);%*model.annotation.defaultLCB;
            model.cofaratios.uratio(Nr_R,1)                = data_raw.data.Metabolites(i,2);%*model.annotation.defaultUCB;
        end
    end
%% Create S-matrices
    for Variant = 1:3
        S = CreateS(model,Variant);
% Add S-matrix to model structure
        if Variant == 1
            model.Variant_Frc.rawS = S;
        elseif Variant == 2
            model.Variant_Suc.rawS = S;
        elseif Variant == 3
            model.Variant_For.rawS = S;
        end
    end
end

function Smatrix = CreateS(model,Variant)
% Initialize the stoichiometric matrix
    
    if Variant == 1
        reaction_string = model.Variant_Frc.reactions.equations;
        rxns            = model.Variant_Frc.reactions.rxns;
    elseif Variant == 2
        reaction_string = model.Variant_Suc.reactions.equations;
        rxns            = model.Variant_Suc.reactions.rxns;
    elseif Variant == 3
        reaction_string = model.Variant_For.reactions.equations;
        rxns            = model.Variant_For.reactions.rxns;
    end
    
    mets            = model.metabolites.mets;
    

    S = zeros(numel(mets),numel(reaction_string));
    reaction_string=strrep(reaction_string,' + ', '¤');
    
% Loop through the reaction_string and add the info to the S matrix
    for i = 1:numel(reaction_string)
    % Start by finding the position of the (=> or <=>)
        arrowIndex = strfind(reaction_string{i},' <=> ');
    
    % Split reactions into reactants and products
        substrates  = regexp(reaction_string{i}(1:arrowIndex-1),'¤','split');
        products    = regexp(reaction_string{i}(arrowIndex+5:end),'¤','split');
    
    % If the splitting character is at the end (if exchange rxns), then an
    % empty string will exist together with the real ones. Remove it
        substrates(cellfun(@isempty,substrates))=[];
        products(cellfun(@isempty,products))=[];

    % A vector where an element is -1 is the corresponding metabolite is a
    % reactant and 1 if it's a product
        multiplyWith=[ones(numel(substrates),1)*-1; ones(numel(products),1)];

        metabolites=[substrates products];

    % Now loop through the metabolites and see if the metabolite has a coefficient 
    % (it will look as 'number name')
            for j=1:numel(metabolites)
                space=strfind(metabolites{j},' ');
        
                if isempty(space)
                %No coefficient
                    coeff=1;
                    name=metabolites{j};
                else
                    coeff=str2double(metabolites{j}(1:space(1)));
            
                %If it was not a coefficiant
                    if isnan(coeff)
                        coeff=1;
                        name=metabolites{j};
                    else
                        name=metabolites{j}(space+1:end);
                    end
                end
        
            %Find the name in the mets list
            %[a b]=ismember(name,mets);
                b=find(strcmp(name,mets),1);
        
                if any(b)
                    S(b,i)=S(b,i)+coeff*multiplyWith(j);
                else
                    if isempty(rxns)
                        dispEM(['Could not find metabolite ' name ' in metabolite list']);
                    else
                        dispEM(['The metabolite "' name '" in reaction ' rxns{i} ' was not found in the metabolite list']);
                    end
                end
            end
    end
    Smatrix = S;
end