%% Get Driving Forces normalized to MDF and Sum of Driving forces around acetyl-CoA

%% Nr 3 in sequence of "after Sampling"


%% Read in Data
addpath('Data/')
load('ThermoCup_DF_Data_All.mat');

DF_Data = DF_Data_Frc;
% DF_Data = DF_Data_Suc;
% DF_Data = DF_Data_For;


%% Prepare Data
max_Nr_EFM = length(DF_Data);
runvar=max_Nr_EFM:-1:1;

for i=1:max_Nr_EFM
    if isempty(DF_Data(runvar(i)).MDF)
        DF_Data(runvar(i))=[];
    end
end

max_Nr_EFM = length(DF_Data);

%% Get XFPK and aclAB Identifier
for i = 1:max_Nr_EFM
    if ~isempty(MJanasch_FindIndex(DF_Data(i).Rxn_Names,'AclAB'))
        DF_Data(i).Strain = 'AclAB';
    elseif ~isempty(MJanasch_FindIndex(DF_Data(i).Rxn_Names,'XFPK1'))
        DF_Data(i).Strain = 'XFPK';
    elseif ~isempty(MJanasch_FindIndex(DF_Data(i).Rxn_Names,'XFPK2'))
        DF_Data(i).Strain = 'XFPK';
    else
        DF_Data(i).Strain = 'WT';
    end
end

%% Remove PhaCE
for i = 1:max_Nr_EFM
    PhaCE_Index = MJanasch_FindIndex(DF_Data(i).Rxn_Names,'PhaCE');
    DF_Data(i).DF(PhaCE_Index) = [];
    %DF_Data(i).MAD(PhaCE_Index) = [];
    DF_Data(i).Rxn_Names(PhaCE_Index) = [];
end


DF_Data_for_Plot = DF_Data;
%% Get which reaction is closest to MDF
for i = 1:max_Nr_EFM
    [DF_Data_for_Plot(i).DF,sorted_Indexes] = sort(DF_Data(i).DF);
    DF_Data_for_Plot(i).Rxn_Names = DF_Data_for_Plot(i).Rxn_Names(sorted_Indexes);
    DF_Data_for_Plot(i).DF_Norm = DF_Data_for_Plot(i).MDF./DF_Data_for_Plot(i).DF;
end

%% Get names of reactions that are used at least once
Rxn_Names = DF_Data_for_Plot(1).Rxn_Names;
for i = 1:max_Nr_EFM-1
    Rxn_Names = [Rxn_Names;DF_Data_for_Plot(i+1).Rxn_Names];
    Rxn_Names = unique(Rxn_Names);
end

DF_for_Plot_Matrix = zeros(length(Rxn_Names),max_Nr_EFM);
DF_for_Plot_Matrix_Norm = zeros(length(Rxn_Names),max_Nr_EFM);

%% Associate each reaction with its DF_Norm for each EFM
for Nr_EFM = 1:max_Nr_EFM
    for Nr_Rxn_Individual = 1:length(DF_Data_for_Plot(Nr_EFM).Rxn_Names)
        Index = MJanasch_FindIndex(Rxn_Names,DF_Data_for_Plot(Nr_EFM).Rxn_Names(Nr_Rxn_Individual));
        if ~isempty(Index)
            DF_for_Plot_Matrix(Index,Nr_EFM) = DF_Data_for_Plot(Nr_EFM).DF(Nr_Rxn_Individual);
            DF_for_Plot_Matrix_Norm(Index,Nr_EFM) = DF_Data_for_Plot(Nr_EFM).DF_Norm(Nr_Rxn_Individual);
        end
    end
    MDF(Nr_EFM)     = DF_Data_for_Plot(Nr_EFM).MDF;
    Strain{Nr_EFM}  = DF_Data_for_Plot(Nr_EFM).Strain;
end

%% Sort after MDF
[sorted_MDF,sorted_Indexes_MDF] = sort(MDF);

DF_for_Plot_Matrix = DF_for_Plot_Matrix(:,sorted_Indexes_MDF);
DF_for_Plot_Matrix_Norm = DF_for_Plot_Matrix_Norm(:,sorted_Indexes_MDF);
sorted_Strain = Strain(sorted_Indexes_MDF);


%% Reactions of Interest
% Reactions producing AcCoA
Reactions_Source = {'PDH';
                    'ACS';
                    'PTA_rev';
                    'AclAB'};

% Reactions consuming AcCoA
Reactions_Sink = {'PhaA';
                    'PTA';
                    'CITS';
                    'MAS'};

%% Get Sums
% For each EFM, get the sum of sources and sums of sinks, and their sum,
% and then display all (see if there is a difference anyway, might be...

Inds_Sources    = ThermoCup_FindIndex(Rxn_Names,Reactions_Source);
Inds_Sinks      = ThermoCup_FindIndex(Rxn_Names,Reactions_Sink);
w = 0;
xf = 0;
a = 0;
for p=1:Nr_EFM
    if strcmp(sorted_Strain{p},'XFPK')
        xf=xf+1;
        Sources_xfpk(xf)  = sum(DF_for_Plot_Matrix(Inds_Sources,p));
        Sinks_xfpk(xf)    = sum(DF_for_Plot_Matrix(Inds_Sinks,p));
        Sum_DF_xfpk(xf)   = Sources_xfpk(xf)-Sinks_xfpk(xf);
    elseif strcmp(sorted_Strain{p},'WT')
        w=w+1;
        Sources_wt(w)  = sum(DF_for_Plot_Matrix(Inds_Sources,p));
        Sinks_wt(w)    = sum(DF_for_Plot_Matrix(Inds_Sinks,p));
        Sum_DF_wt(w)   = Sources_wt(w)-Sinks_wt(w);
    else
        a=a+1;
        Sources_aclAB(a)  = sum(DF_for_Plot_Matrix(Inds_Sources,p));
        Sinks_aclAB(a)    = sum(DF_for_Plot_Matrix(Inds_Sinks,p));
        Sum_DF_aclAB(a)   = Sources_aclAB(a)-Sinks_aclAB(a);
    end
end

% [sorted_Sum_DF,i_sorted_Sum_DF] = sort(Sum_DF);
% sorted_Strain_DF = sorted_Strain(i_sorted_Sum_DF);