%% Read in Range data and transform for R output


%% Nr 2 in sequence of "after Sampling"


%% Read in model

addpath('Model/')
modelFile = '/Model/Cupriavidus_Core_XII_PHB.xlsx';
CupriaCore = ThermoCup_ReadExcel(modelFile);


%% Choose Substrate
Substrate = 'Frc';
%Substrate = 'Suc';
%Substrate = 'For';




model_Frc = CupriaCore.Variant_Frc;

Rxn_Names = model_Frc.reactions.rxns;
Met_Names = CupriaCore.metabolites.mets;

Nr_Rxns = length(Rxn_Names);
Nr_Mets = length(Met_Names);

%% Add path to Metabolite range files
addpath('Output_for_Sampling_Frc/')


%% Modify Met_Names

for i = 1:Nr_Mets
    Met_Names{i} = ['M_',Met_Names{i}];
end



%% Read in Data

a = dir('Output_for_Sampling_Frc/MDF/*.txt');
Nr_EFMs = size(a,1);
disp(['Number of analyzed EFMs: ',num2str(Nr_EFMs)])

Min_Value = zeros(Nr_Mets,Nr_EFMs);
Max_Value = zeros(Nr_Mets,Nr_EFMs);
MDF_ALL = zeros(Nr_EFMs,1);
tic
for EFMs_Index = 1:Nr_EFMs % Loop through EFMs
%     Medians_Vector  = zeros(Nr_Rxns_Pruned,1);
%     MADs_Vector     = zeros(Nr_Rxns_Pruned,1);
    
        clear Names_EFMs
        clear Met_Range_Value
        clear Met_Ranges_EFMs
        clear MDF_EFMs
        clear Met_Name_EFM

        % Load Names
        Name_file_name = ['Output_for_Sampling_Frc/Name_References/Name_References_EFM_Nr_',int2str(EFMs_Index),'_Frc.txt'];
        Names_EFMs = importdata(Name_file_name);
        r=0;
        for q = 1:length(Names_EFMs)
            if strncmp(Names_EFMs(q),'M_',2)
                r = r+1;
                Met_Name_EFM(r) = erase(Names_EFMs(q),'M_');
            end
        end



        % Load MDF
        MDF_file_name = ['Output_for_Sampling_Frc/MDF/MDF_EFM_Nr_',int2str(EFMs_Index),'_Frc.txt'];
        MDF_EFMs = importdata(MDF_file_name);

        MDF_ALL(EFMs_Index,1) = MDF_EFMs;

        % Load Met_Range file
        Met_Range_File_Name = ['Output_for_Sampling_Frc/Met_Ranges/MetRange_EFM_Nr_',int2str(EFMs_Index),'_Frc.txt'];
        Met_Ranges_EFMs = importdata(Met_Range_File_Name);

        for p = 1:length(Met_Ranges_EFMs)
            Met_Range_Min(p,1) = str2double(extractBetween(Met_Ranges_EFMs(p),'[',' '));
            Met_Range_Max(p,2) = str2double(extractBetween(Met_Ranges_EFMs(p),' ',']'));
        end


        Met_Range_Result(EFMs_Index).MDF = MDF_EFMs;
        Met_Range_Result(EFMs_Index).Met_Range.Min = Met_Range_Min;
        Met_Range_Result(EFMs_Index).Met_Range.Max = Met_Range_Max;    
        Met_Range_Result(EFMs_Index).Met_Names = Met_Name_EFM;

        if rem(EFMs_Index,1000) == 0
            disp(['Number of EFMs analyzed: ',num2str(EFMs_Index)])
            toc
        end
end

MDF = MDF_ALL;



%% Read in Data
addpath('Data/')
load('ThermoCup_All_Enumerated_EFMs.mat');

mnet = mnet_Frc;


%% Remove empty rows in vector
% max_Nr_EFM = length(Met_Range_Result);
% runvar=max_Nr_EFM:-1:1;
% 
% for i=1:max_Nr_EFM
%     if isempty(Met_Range_Result(runvar(i)).MDF)
%         Met_Range_Result(runvar(i))=[];
%     end
% end
% 
% max_Nr_EFM = length(Met_Range_Result);
% 
% 
% for i=1:max_Nr_EFM
%     MDF(i) = Met_Range_Result(i).MDF;
% end
% 
% sorted_MDF = sort(MDF);
% x=1:1:max_Nr_EFM;


%% loop through files and check if xfpk, aclAB or WT



%% Put met ranges into right format
clear Met_Range
Met_Range.Min = [];
Met_Range.Max = [];

MDF_yes = 0;

for p = 1:length(MDF) % Loop through each EFM
    if MDF(p) >0
        MDF_yes = MDF_yes + 1;
        clear Met_List
        Met_List = string(transpose(Met_Range_Result(p).Met_Names)); % Get metabolites of EFM
        %Met_List(Met_List=="PHB[c]") =[]; % Remove PHB from metabolite list, sometimes in there double? Not supposed to change anyhow                
        for u = 1:length(mnet.metaboliteNames) % loop through all metabolites in the model
            Ind = MJanasch_FindIndex(Met_List,mnet.metaboliteNames{u}); % Find index of metabolite in EFM_metabolite list
            if ~isempty(Ind) % If metabolite used in EFM
%                 Met_Range(u).Min(MDF_yes)= EVA(MDF_yes).Results(Ind,1).Conc_Out;
%                 Met_Range(u).Max(MDF_yes)= EVA(MDF_yes).Results(Ind,2).Conc_Out;
                Met_Range(u).Min(MDF_yes)= Met_Range_Result(p).Met_Range.Min(Ind);
                Met_Range(u).Max(MDF_yes)= Met_Range_Result(p).Met_Range.Max(Ind,2);
                
                
            else
                Met_Range(u).Min(MDF_yes)= 1000;
                Met_Range(u).Max(MDF_yes)= 1000;
            end       
        end
    end
end

% Remove all 0 entries
for u = 1:length(mnet.metaboliteNames)
    Indexes=MJanasch_FindIndex(Met_Range(u).Max,0);
    if ~isempty(Indexes)
        Met_Range(u).Min(Indexes) = [];
        Met_Range(u).Max(Indexes) = [];
    end
%     Met_Range(u).MeanMin = mean(Met_Range(u).Min);
%     Met_Range(u).MeanMax = mean(Met_Range(u).Max);
%      
end



num_EFM_Frc = length(MDF);
x_Frc=1:1:num_EFM_Frc;
Nr_Mets = length(Met_Range);


%% Sort Met Ranges after MDF and DF
load('Indices_DoubleSorted_Frc.mat')

MetStrats = unique(MDF);
Nr_MetStrats = length(MetStrats);

MetStrat_EFM_Count = 0;
indices_MetStrat_all = [];

[MDF_sorted_Frc,i_Frc_MDF] = sort(MDF);



% MetStrat_EFM_Count_Additive = 0;
% for i = 1:Nr_MetStrats              % For each metabolic strategies
%     
%     % Find All EFMs of that strategy
%     MetStrat_Indices = ThermoCup_FindIndex(MDF_sorted_Frc,MetStrats(i));
%     Nr_EFMs_in_MetStrat = length(MetStrat_Indices);
%     
% 
%     for j = 1:Nr_EFMs_in_MetStrat
%         Indices_DoubleSorted_Frc_MOD(j+MetStrat_EFM_Count,1) = Indices_DoubleSorted_Frc(j+MetStrat_EFM_Count)+MetStrat_EFM_Count_Additive(i);
%     end
%     
% MetStrat_EFM_Count = MetStrat_EFM_Count+Nr_EFMs_in_MetStrat;
% MetStrat_EFM_Count_All(i)=Nr_EFMs_in_MetStrat;
% MetStrat_EFM_Count_Additive(i+1) = MetStrat_EFM_Count_Additive(i) + MetStrat_EFM_Count_All(i);
% end




clear Met_Range_sorted_MDF

%% Put Metabolite Ranges into sorted order
for i = 1:Nr_Mets
    Met_Range(i).Min = Met_Range(i).Min(i_Frc_MDF);
    Met_Range(i).Max = Met_Range(i).Max(i_Frc_MDF);
    
    Met_Range_sorted_MDF(i).Min = log10(Met_Range(i).Min(Indices_DoubleSorted_Frc));
    Met_Range_sorted_MDF(i).Max = log10(Met_Range(i).Max(Indices_DoubleSorted_Frc));
end

Met_Range_sorted_MDF_All_Frc = Met_Range_sorted_MDF;
%% Prepare for R

x_Frc=1:1:length(MDF_sorted_Frc);
To_plot_Frc_Min(:,1) = x_Frc;
To_plot_Frc_Max(:,1) = x_Frc;
To_plot_Frc_Min(:,2) = MDF_sorted_Frc;
To_plot_Frc_Max(:,2) = MDF_sorted_Frc;


%% Put metabolite ranges
for i=1:Nr_Mets
%     To_plot_Frc_Min(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Min);
%     To_plot_Frc_Max(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Max);
%     
%     To_plot_Frc_Min(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Min);
%     To_plot_Frc_Max(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Max);
    
    To_plot_Frc_Min(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Min);
    To_plot_Frc_Max(:,i+2) = transpose(Met_Range_sorted_MDF_All_Frc(i).Max);

    %% Make Header
    Met_Names_Mod{i} = replace(CupriaCore.metabolites.mets{i},'[','_');
    Met_Names_Mod{i} = erase(Met_Names_Mod{i},']');
    
end    





Header = ['EFM','MDF',Met_Names_Mod];
% To_plot_Frc_Min(2:end+1,:) = To_plot_Frc_Min;
% To_plot_Frc_Min(1,:) = Header;
%% Save data

save('To_plot_Frc_Min.txt','To_plot_Frc_Min', '-ascii');
save('To_plot_Frc_Max.txt','To_plot_Frc_Max', '-ascii');

% save('To_plot_Suc_Min.txt','To_plot_Suc_Min', '-ascii');
% save('To_plot_Suc_Max.txt','To_plot_Suc_Max', '-ascii');

% save('To_plot_For_Min.txt','To_plot_For_Min', '-ascii');
% save('To_plot_For_Max.txt','To_plot_For_Max', '-ascii');