%% Checks percentage reaction usage over all EFMs

%%%------- Requires to run ThermoCup_EnumerateEFMs.m first to get filtered EFM counts!!!
% Depending on job, ThermoCup_EnumerateEFMs.m needs to be run finding only
% EFMs producing biomass, PHB or both

clear EFM_Count_allpos
clear EFM_Count_Sum
clear EFM_Count_Split
clear Reaction_Name_Rev
clear EFM_Count_total

%% Load Data framework
addpath('Data/')
load('ThermoCup_II_All_Enumerated_EFMs.mat');

mnet = mnet_Frc;
% mnet = mnet_For;
Nr_Rxns = length(mnet.reactionNames);

%% Load Data
% Load the dataset for the corresponding substrate, created from ThermoCup_EnumerateEFMs.m
 load('ThermoCup_II_Yields_MDF_Frc_Bio.mat');
% load('ThermoCup_II_Yields_MDF_Frc_PHB.mat');

% load('ThermoCup_II_Yields_MDF_For_Bio.mat');
%load('ThermoCup_II_Yields_MDF_For_PHB.mat');


%% Split reversible reactions
r=0;

% Decide for which strain's reaction usage is to be analyzed
EFM_Count_total = EFM_Count;
% EFM_Count_total = EFM_Count_xfpk;
% EFM_Count_total = EFM_Count_aclAB;

for i = 1:Nr_Rxns                                       % Loop through reactions
    r=r+1;                                              % Reaction index
    Reaction_Name_Rev(r) = mnet.reactionNames(i);       % Forward reaction name
    if any(EFM_Count_total(i,:)<0)                      % If reaction is used in reverse
        EFM_Count_Split(r,:) = EFM_Count_total(i,:);    % All directionalities first copied
            for j = 1:size(EFM_Count_total,2)           % loop through EFMs
                if EFM_Count_Split(r,j) <0
                    EFM_Count_Split(r,j) = 0;
                end
            end
        r=r+1;
        EFM_Count_Split(r,:) = -EFM_Count_total(i,:);
            for j = 1:size(EFM_Count_total,2)
                if EFM_Count_Split(r,j) < 0
                    EFM_Count_Split(r,j) = 0;
                end
            end
        Reaction_Name_Rev(r) = {[mnet.reactionNames{i},'-rev']};
    else
        EFM_Count_Split(r,:) = EFM_Count_total(i,:);
    end
end

%% Get percentage of reaction usage
numtot = size(EFM_Count_total,2);
EFM_Count_allpos = EFM_Count_Split;
Reaction_Name_Rev = transpose(Reaction_Name_Rev);
Nr_Rxns_Rev = length(Reaction_Name_Rev);

for ii = 1:numtot
    for xy = 1:Nr_Rxns_Rev
        if EFM_Count_Split(xy,ii) <0
            EFM_Count_allpos(xy,ii) = 1;
        end
    end
end

for xyz = 1:Nr_Rxns_Rev
    EFM_Count_Sum(xyz,1) = 100*sum(EFM_Count_allpos(xyz,:))/numtot;
end