
function model = ThermoCup_ReadExcel_FBA(fileName)


% fileName = '/Users/markus.janasch/Dropbox/Cyanogroup/Project Documents/Thermodynamics/Cupriavidus_Markus/AcCoA_MVA_Project/Input/Cupriavidus_Core.xlsx';
data_raw = importdata(fileName);
%     model.id            = data_raw.textdata.Info{2,1};
%     model.description   = data_raw.textdata.Info{2,2};

% Reactions
    model.rxns        = data_raw.textdata.Reactions(2:end,1);
    model.rxnNames    = data_raw.textdata.Reactions(2:end,3);
    model.equations   = data_raw.textdata.Reactions(2:end,4);
    model.lb          = data_raw.data.Reactions(:,3);
    model.ub          = data_raw.data.Reactions(:,4);
    model.eccodes     = data_raw.textdata.Reactions(2:end,6);
    model.dG0         = data_raw.data.Reactions(:,1);
    
    for i = 1:length(model.rxns)
        if model.lb(i) < 0
            model.rev(i,1)  = 1;
        else
            model.rev(i,1)  = 0;
        end
    end
    model.c = zeros(length(model.rxns),1);
    
%% Metabolites
    model.mets      = data_raw.textdata.Metabolites(2:end,1);
    model.metNames  = data_raw.textdata.Metabolites(2:end,2);
    model.lconc     = data_raw.data.Metabolites(:,1);%*model.annotation.defaultLCB;
    model.uconc     = data_raw.data.Metabolites(:,2);%*model.annotation.defaultUCB;
    model.b         = zeros(length(model.mets),1);
    
%% Create S-matrix
    

    reaction_string = model.equations;
    rxns            = model.rxns;
    mets            = model.mets;
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
    model.S = Smatrix;
end