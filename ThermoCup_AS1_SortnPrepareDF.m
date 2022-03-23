%% Sort EFMs after MDF and inside each MDF-value

%% Nr 1 in sequence of "after Sampling"


%Substrate = 'Frc';
%Substrate = 'Suc';
Substrate = 'For';


a = dir(['Output_for_Sampling_',Substrate,'/',Substrate,'_PHB_Median/*.txt']);
Nr_EFMs = size(a,1);
disp(['Number of analyzed EFMs: ',num2str(Nr_EFMs)])

tic
for EFMs_Index = 1:Nr_EFMs % Loop through EFMs

        clear Names_EFMs
        clear DF_Medians
        clear DF_MADs
        clear MDF_EFMs
        clear Rxn_Name_EFM

%% Load Names
        Name_file_name = ['Output_for_Sampling_',Substrate,'/Name_References/Name_References_EFM_Nr_',int2str(EFMs_Index),'_',Substrate,'.txt'];
        Names_EFMs = importdata(Name_file_name);
        r=0;
        for q = 1:length(Names_EFMs)
            if strncmp(Names_EFMs(q),'R_',2)
                r = r+1;
                Rxn_Name_EFM(r,1) = erase(Names_EFMs(q),'R_');
            end
        end

%% Load MDFs
        MDF_file_name = ['Output_for_Sampling_',Substrate,'/MDF/MDF_EFM_Nr_',int2str(EFMs_Index),'_',Substrate,'.txt'];
        MDF_EFMs = importdata(MDF_file_name);

        %MDF_ALL(EFMs_Index) = MDF_EFMs;

%% Load DFs
        DF_Median_File_Name = ['Output_for_Sampling_',Substrate,'/',Substrate,'_PHB_Median/Sampling_Results_MEDIAN_EFM_Nr_',int2str(EFMs_Index),'_PHB_',Substrate,'.txt'];
        DF_Median_EFMs = importdata(DF_Median_File_Name);

%% Load MAD DFs
        DF_MAD_File_Name = ['Output_for_Sampling_',Substrate,'/',Substrate,'_PHB_MAD/Sampling_Results_MAD_EFM_Nr_',int2str(EFMs_Index),'_PHB_',Substrate,'.txt'];
        DF_MAD_EFMs = importdata(DF_MAD_File_Name);

        
%% Combine data in one structure

        
        MDF_Data(EFMs_Index) = MDF_EFMs;
%         Rxn_Names_Data(EFMs_Index) = Rxn_Name_EFM;
%         DF_Median_Data(EFMs_Index) = DF_Median_EFMs;
%         DF_MAD_Data(EFMs_Index) = DF_MAD_EFMs;
        
        DF_Data(EFMs_Index).MDF = MDF_EFMs;
        DF_Data(EFMs_Index).Rxn_Names = Rxn_Name_EFM;
        DF_Data(EFMs_Index).DF = DF_Median_EFMs;    
        DF_Data(EFMs_Index).MAD = DF_MAD_EFMs;

        if rem(EFMs_Index,1000) == 0
            disp(['Number of EFMs analyzed: ',num2str(EFMs_Index)])
            toc
        end
end

[sorted_MDF,i_sorted] = sort(MDF_Data);


for i = 1:Nr_EFMs
    DF_Data_new(i).MDF = MDF_Data(i_sorted(i));
    DF_Data_new(i).Rxn_Names = DF_Data(i_sorted(i)).Rxn_Names;
    DF_Data_new(i).DF = DF_Data(i_sorted(i)).DF;
    DF_Data_new(i).MeanDF = mean(DF_Data_new(i).DF);
    MeanDF(i,1) = mean(DF_Data_new(i).DF);
    DF_Data_new(i).MAD = DF_Data(i_sorted(i)).MAD;
end

%% Find Metabolic Strategies (i.e. all different MDFs)

MetStrats = unique(sorted_MDF);
Nr_MetStrats = length(MetStrats);

MetStrat_Count = 0;
indices_MetStrat_all = [];
for i = 1:Nr_MetStrats              % For each metabolic strategies
    
    % Find All EFMs of that strategy
    MetStrat_Indices = ThermoCup_FindIndex(sorted_MDF,MetStrats(i));
    Nr_EFMs_in_MetStrat = length(MetStrat_Indices);
    Mean_DF_in_MetStrat = MeanDF(MetStrat_Indices);
    % Sort EFMs by Mean
    [sorted_Mean_DF,ind_sorted_MeanDF] = sort(Mean_DF_in_MetStrat);
    sorted_MetStrat_Indices = MetStrat_Indices(ind_sorted_MeanDF);
    
    indices_MetStrat_all = [indices_MetStrat_all;ind_sorted_MeanDF+MetStrat_Count];
    
    for j = 1:Nr_EFMs_in_MetStrat
        DF_Data_Final(j+MetStrat_Count).MDF        = MetStrats(i);
        DF_Data_Final(j+MetStrat_Count).Rxn_Names  = DF_Data_new(sorted_MetStrat_Indices(j)).Rxn_Names;
        DF_Data_Final(j+MetStrat_Count).DF         = DF_Data_new(sorted_MetStrat_Indices(j)).DF;
        DF_Data_Final(j+MetStrat_Count).MeanDF     = DF_Data_new(sorted_MetStrat_Indices(j)).MeanDF;
        DF_Data_Final(j+MetStrat_Count).MAD        = DF_Data_new(sorted_MetStrat_Indices(j)).MAD;
    end
    
    
MetStrat_Count = MetStrat_Count+Nr_EFMs_in_MetStrat;
end
eval(['DF_Data_',Substrate,' = DF_Data_Final;'])
eval(['Indices_DoubleSorted_',Substrate,' = indices_MetStrat_all;'])

cd(['Output_for_Sampling_',Substrate,'/']);
save(['Sorted_EFMs_',Substrate,'_DF.mat'],['DF_Data_',Substrate])
save(['Indices_DoubleSorted_',Substrate,'.mat'],['Indices_DoubleSorted_',Substrate])