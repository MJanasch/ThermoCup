%% DF for formate map

%% Nr 4 in sequence of "after Sampling"


load('Output_for_Sampling_For/Sorted_EFMs_For_DF.mat')

addpath('Model/')
modelFile = '/Model/Cupriavidus_Core_XII_PHB.xlsx';
CupriaCore = ThermoCup_ReadExcel(modelFile);


%% Find Metabolic Strategies (i.e. all different MDFs)
for i=1:length(DF_Data_For)
    sorted_MDF(i,1) = DF_Data_For(i).MDF;
end

MetStrats = unique(sorted_MDF);
Nr_MetStrats = length(MetStrats);
MetStrat_Count = 0;



Rxn_Names = DF_Data_For(1).Rxn_Names;
for i = 1:length(DF_Data_For)-1
    Rxn_Names = [Rxn_Names;DF_Data_For(i+1).Rxn_Names];
    Rxn_Names = unique(Rxn_Names);
end




for i = 1:Nr_MetStrats              % For each metabolic strategies
    
    % Find All EFMs of that strategy
    MetStrat_Indices = ThermoCup_FindIndex(sorted_MDF,MetStrats(i));
    Nr_EFMs_in_MetStrat = length(MetStrat_Indices);
    
    
    
    
    
    % Loop through all EFMs of a Metabolic Strategy and collect all DF
    for j=1:Nr_EFMs_in_MetStrat
        for k=1:length(DF_Data_For(MetStrat_Indices(j)).Rxn_Names)
            Rxn_Index = MJanasch_FindIndex(Rxn_Names,DF_Data_For(MetStrat_Indices(j)).Rxn_Names(k));
            if ~isempty(Rxn_Index)
                DF_Reaction(Rxn_Index,j) = DF_Data_For(MetStrat_Indices(j)).DF(k);
            else
                DF_Reaction(Rxn_Index,j) = 0;
            end
        end
            
    end
    

MetStrat_Count = MetStrat_Count+Nr_EFMs_in_MetStrat;
eval(['DF_Reaction_MetStrat_',int2str(i),' = DF_Reaction;']);
clear DF_Reaction
end

for p=1:55
    if sum(DF_Reaction_MetStrat_1(p,:)) ~=0
        Stat_DF_1(p,1) = mean(nonzeros(DF_Reaction_MetStrat_1(p,:)));
        Stat_DF_1(p,2) = std(nonzeros(DF_Reaction_MetStrat_1(p,:)));
        Stat_DF_1(p,3) = median(nonzeros(DF_Reaction_MetStrat_1(p,:)));
        Stat_DF_1(p,4) = mad(nonzeros(DF_Reaction_MetStrat_1(p,:)),1);
    else
        Stat_DF_1(p,1) = 0;
        Stat_DF_1(p,2) = 0;
        Stat_DF_1(p,3) = 0;
        Stat_DF_1(p,4) = 0;
    end
end

for p=1:55
    if sum(DF_Reaction_MetStrat_2(p,:)) ~=0
        Stat_DF_2(p,1) = mean(nonzeros(DF_Reaction_MetStrat_2(p,:)));
        Stat_DF_2(p,2) = std(nonzeros(DF_Reaction_MetStrat_2(p,:)));
        Stat_DF_2(p,3) = median(nonzeros(DF_Reaction_MetStrat_2(p,:)));
        Stat_DF_2(p,4) = mad(nonzeros(DF_Reaction_MetStrat_2(p,:)),1);
    else
        Stat_DF_2(p,1) = 0;
        Stat_DF_2(p,2) = 0;
        Stat_DF_2(p,3) = 0;
        Stat_DF_2(p,4) = 0;
    end
end
