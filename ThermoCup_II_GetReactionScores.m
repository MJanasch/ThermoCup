%% Calculate Yield- and MDF-Scores

%%%------- Requires to run ThermoCup_EnumerateEFMs.m first to get filtered EFM counts!!!
% Depending on job, ThermoCup_EnumerateEFMs.m needs to be run finding only
% EFMs producing biomass, PHB or both, and name the EFM_Count data and
% yield data accordingly


%% Read in Data
addpath('Data/')
load('ThermoCup_II_All_Enumerated_EFMs.mat');

% Load the dataset for the corresponding substrate, created from ThermoCup_EnumerateEFMs.m
load('ThermoCup_II_Yields_MDF_Frc.mat');
% load('ThermoCup_II_Yields_MDF_For.mat');


%% Choose a substrate
% Substrate = 'Frc';
Substrate = 'For';

if exist('Substrate','var') == 1
    if strcmp(Substrate,'Frc')
        mnet = mnet_Frc;
        EFM_Count_total     = [EFM_Count EFM_Count_aclAB EFM_Count_xfpk];
        yield_Biomass_total_All = [yield_Biomass yield_Biomass_aclAB yield_Biomass_xfpk];
        yield_PHB_total_All     = [yield_PHB yield_PHB_aclAB yield_PHB_xfpk];
        MDF_total = [MDF_yes MDF_aclAB_yes MDF_xfpk_yes];
        disp('Substrate is Frc')
    elseif strcmp(Substrate,'For')
        mnet = mnet_For;
        EFM_Count_total     = [EFM_Count EFM_Count_xfpk];
        yield_Biomass_total_All = [yield_Biomass yield_Biomass_xfpk];
        yield_PHB_total_All     = [yield_PHB yield_PHB_xfpk];
        MDF_total = [MDF_yes MDF_xfpk_yes];
        disp('Substrate is For')
    else
        disp('Choose a right substrate')
    end
else
    disp('Define a substrate, please choose between "Frc" or "For".')
    return
end


Nr_Rxns = length(mnet.reactionNames);


%% Get right reaction usage

i_bio = 0;
i_PHB = 0;

for i=1:length(EFM_Count_total)
    if EFM_Count_total(4,i) ~=0
        i_bio = i_bio + 1;
        EFM_Count_total_Biomass(:,i_bio) = EFM_Count_total(:,i);
        yield_Biomass_total(i_bio) = yield_Biomass_total_All(i);
    end
end


for i=1:length(EFM_Count_total)
    if EFM_Count_total(9,i) ~=0
        i_PHB = i_PHB + 1;
        EFM_Count_total_PHB(:,i_PHB) = EFM_Count_total(:,i);
        yield_PHB_total(i_PHB) = yield_PHB_total_All(i);
    end
end

EFM_Count_total_MDF = EFM_Count_total;



%% Split reversible reactions

% Biomass yield
r=0;
for i = 1:Nr_Rxns
    r=r+1;
    Reaction_Name_Biomass(r) = mnet.reactionNames(i);
    if any(EFM_Count_total_Biomass(i,:)<0)
        EFM_Count_Split_Biomass(r,:) = EFM_Count_total_Biomass(i,:);
            for j = 1:length(yield_Biomass_total)
                if EFM_Count_Split_Biomass(r,j) <0
                    EFM_Count_Split_Biomass(r,j) = 0;
                end
            end
        r=r+1;
        EFM_Count_Split_Biomass(r,:) = -EFM_Count_total_Biomass(i,:);
            for j = 1:length(yield_Biomass_total)
                if EFM_Count_Split_Biomass(r,j) < 0
                    EFM_Count_Split_Biomass(r,j) = 0;
                end
            end
        Reaction_Name_Biomass(r) = {[mnet.reactionNames{i},'-rev']};
    else
        EFM_Count_Split_Biomass(r,:) = EFM_Count_total_Biomass(i,:);
    end
end
Reaction_Name_Before_Biomass = transpose(Reaction_Name_Biomass);


count = length(Reaction_Name_Before_Biomass):-1:1;
for i = 1:length(Reaction_Name_Before_Biomass)
    if sum(EFM_Count_Split_Biomass(count(i),:)) == 0
        EFM_Count_Split_Biomass(count(i),:) = [];
        Reaction_Name_Before_Biomass(count(i)) = [];
    end
end



% PHB yield
r=0;
for i = 1:Nr_Rxns
    r=r+1;
    Reaction_Name_PHB(r) = mnet.reactionNames(i);
    if any(EFM_Count_total_PHB(i,:)<0)
        EFM_Count_Split_PHB(r,:) = EFM_Count_total_PHB(i,:);
            for j = 1:length(yield_PHB_total)
                if EFM_Count_Split_PHB(r,j) <0
                    EFM_Count_Split_PHB(r,j) = 0;
                end
            end
        r=r+1;
        EFM_Count_Split_PHB(r,:) = -EFM_Count_total_PHB(i,:);
            for j = 1:length(yield_PHB_total)
                if EFM_Count_Split_PHB(r,j) < 0
                    EFM_Count_Split_PHB(r,j) = 0;
                end
            end
        Reaction_Name_PHB(r) = {[mnet.reactionNames{i},'-rev']};
    else
        EFM_Count_Split_PHB(r,:) = EFM_Count_total_PHB(i,:);
    end
end
Reaction_Name_Before_PHB = transpose(Reaction_Name_PHB);

count = length(Reaction_Name_Before_PHB):-1:1;
for i = 1:length(Reaction_Name_Before_PHB)
    if sum(EFM_Count_Split_PHB(count(i),:)) == 0
        EFM_Count_Split_PHB(count(i),:) = [];
        Reaction_Name_Before_PHB(count(i)) = [];
    end
end


% MDF
r=0;
for i = 1:Nr_Rxns
    r=r+1;
    Reaction_Name_MDF(r) = mnet.reactionNames(i);
    if any(EFM_Count_total_MDF(i,:)<0)
        EFM_Count_Split_MDF(r,:) = EFM_Count_total_MDF(i,:);
            for j = 1:length(MDF_total)
                if EFM_Count_Split_MDF(r,j) <0
                    EFM_Count_Split_MDF(r,j) = 0;
                end
            end
        r=r+1;
        EFM_Count_Split_MDF(r,:) = -EFM_Count_total_MDF(i,:);
            for j = 1:length(MDF_total)
                if EFM_Count_Split_MDF(r,j) < 0
                    EFM_Count_Split_MDF(r,j) = 0;
                end
            end
        Reaction_Name_MDF(r) = {[mnet.reactionNames{i},'-rev']};
    else
        EFM_Count_Split_MDF(r,:) = EFM_Count_total_MDF(i,:);
    end
end
Reaction_Name_Before_MDF = transpose(Reaction_Name_MDF);

count = length(Reaction_Name_Before_MDF):-1:1;
for i = 1:length(Reaction_Name_Before_MDF)
    if sum(EFM_Count_Split_MDF(count(i),:)) == 0
        EFM_Count_Split_MDF(count(i),:) = [];
        Reaction_Name_Before_MDF(count(i)) = [];
    end
end



%% Calculate Yield and MDF Scores 2.0

total_mean_yield_Biomass = mean(yield_Biomass_total);

total_mean_yield_PHB = mean(yield_PHB_total);

total_mean_MDF = mean(MDF_total);


for p = 1:length(Reaction_Name_Before_Biomass)
    
    Nr_Used_Reactions_Biomass = sum(EFM_Count_Split_Biomass(p,:));
    mean_yield_Biomass(p) = (EFM_Count_Split_Biomass(p,:) * transpose(yield_Biomass_total))/Nr_Used_Reactions_Biomass;
    Biomass_Yield_Score(p) = (round((mean_yield_Biomass(p)/total_mean_yield_Biomass)-1,5))*100;
end

for q = 1:length(Reaction_Name_Before_PHB)
    Nr_Used_Reactions_PHB = sum(EFM_Count_Split_PHB(q,:));
    mean_yield_PHB(q) = (EFM_Count_Split_PHB(q,:) * transpose(yield_PHB_total))/Nr_Used_Reactions_PHB;
    PHB_Yield_Score(q) = (round((mean_yield_PHB(q)/total_mean_yield_PHB)-1,5))*100;
end

for s = 1:length(Reaction_Name_Before_MDF)
    Nr_Used_Reactions_MDF = sum(EFM_Count_Split_MDF(s,:));
    mean_MDF(s) = (EFM_Count_Split_MDF(s,:) * transpose(MDF_total))/Nr_Used_Reactions_MDF;
    MDF_Score(s) = (round((mean_MDF(s)/total_mean_MDF)-1,5))*100;
end


% %% Remove unused reactions
% runvar_back = length(Reaction_Name_Before_Biomass):-1:1;
% for i = 1:length(runvar_back)
%     if Biomass_Yield_Score(runvar_back(i)) == 0 && PHB_Yield_Score(runvar_back(i)) == 0
%         Reaction_Name_Before_Biomass(runvar_back(i))    = [];
%         Reaction_Name_Before_PHB(runvar_back(i))        = [];
%         MDF_Score(runvar_back(i))                       = [];
%         Biomass_Yield_Score(runvar_back(i))             = [];
%         PHB_Yield_Score(runvar_back(i))                 = [];
%     elseif strncmp(Reaction_Name_Before_Biomass(runvar_back(i)),'Ex_H2O',6)
%         Reaction_Name_Before_Biomass(runvar_back(i))    = [];
%         Reaction_Name_Before_PHB(runvar_back(i))        = [];
%         MDF_Score(runvar_back(i))                       = [];
%         Biomass_Yield_Score(runvar_back(i))             = [];
%         PHB_Yield_Score(runvar_back(i))                 = [];
%     elseif strncmp(Reaction_Name_Before_Biomass(runvar_back(i)),'Trp_H2O',6)
%         Reaction_Name_Before_Biomass(runvar_back(i))    = [];
%         Reaction_Name_Before_PHB(runvar_back(i))        = [];
%         MDF_Score(runvar_back(i))                       = [];
%         Biomass_Yield_Score(runvar_back(i))             = [];
%         PHB_Yield_Score(runvar_back(i))                 = [];
%     end
% end

MDF_Score_trans = transpose(MDF_Score);
Biomass_Yield_Score_trans = transpose(Biomass_Yield_Score);
PHB_Yield_Score_trans = transpose(PHB_Yield_Score);
Reaction_Name_trans = Reaction_Name_Before_Biomass;

%dx = 0.5; dy = 0.5;

%% Plot
% figure(1)
% %plot(MDF_Score_trans,Yield_Score_trans,'*k')
% scatter(Biomass_Yield_Score_trans,PHB_Yield_Score_trans,50,MDF_Score_trans,'filled')
% cb = colorbar;
% text(Biomass_Yield_Score_trans+dx,PHB_Yield_Score_trans+dy,Reaction_Name_trans)
% %xlim([-50 30])
% % ylim([-200 500])
% grid on

%% Do Wilcoxon rank sum test and FDR

% Biomass
[All_Directions,Nr_EFMs] = size(EFM_Count_Split_Biomass);
used = 0;
for l=1:All_Directions  % Loop through all reactions (including reverse)
    clear Test_Reaction
    clear Test_Reaction_Null
    clear Test_Reaction_Yes
    if abs(sum(EFM_Count_Split_Biomass(l,:))) ~= Nr_EFMs && ~sum(EFM_Count_Split_Biomass(l,:)) == 0 % if the reaction is not used in all EFMs and is used (e.g. it's not zero)
        used = used + 1;
        Test_Reaction = EFM_Count_Split_Biomass(l,:); % Take the reaction to be tested
        Reac_Bio{used,1} = Reaction_Name_Before_Biomass(l);
        Score_Bio(used,1) = Biomass_Yield_Score_trans(l);
        for q = 1:length(Test_Reaction) 
            if  Test_Reaction(q) == 1
                Test_Reaction_Null(q) = 0;
                Test_Reaction_Yes(q) = yield_Biomass_total(q);
            else
                Test_Reaction_Null(q) = yield_Biomass_total(q);
                Test_Reaction_Yes(q) = 0;
            end
        end
        Test_Reaction_Null(Test_Reaction_Null==0) = [];
        Test_Reaction_Yes(Test_Reaction_Yes==0) = [];
        [pvalue_Bio(used,1),h_Bio(used)] = ranksum(Test_Reaction_Null,Test_Reaction_Yes);
        
    end
        
end
FDR_Biomass = mafdr(pvalue_Bio,'BHFDR','true');


% PHB
[All_Directions,Nr_EFMs] = size(EFM_Count_Split_PHB);
used = 0;
for l=1:All_Directions  % Loop through all reactions (including reverse)
    clear Test_Reaction
    clear Test_Reaction_Null
    clear Test_Reaction_Yes
    if abs(sum(EFM_Count_Split_PHB(l,:))) ~= Nr_EFMs && ~sum(EFM_Count_Split_PHB(l,:)) == 0 % if the reaction is not used in all EFMs and is used (e.g. it's not zero)
        used = used + 1;
        Test_Reaction = EFM_Count_Split_PHB(l,:); % Take the reaction to be tested
        Reac_PHB{used,1} = Reaction_Name_Before_PHB(l);
        Score_PHB(used,1) = PHB_Yield_Score_trans(l);
        for q = 1:length(Test_Reaction) 
            if  Test_Reaction(q) == 1
                Test_Reaction_Null(q) = 0;
                Test_Reaction_Yes(q) = yield_PHB_total(q);
            else
                Test_Reaction_Null(q) = yield_PHB_total(q);
                Test_Reaction_Yes(q) = 0;
            end
        end
        Test_Reaction_Null(Test_Reaction_Null==0) = [];
        Test_Reaction_Yes(Test_Reaction_Yes==0) = [];
        [pvalue_PHB(used,1),h_PHB(used)] = ranksum(Test_Reaction_Null,Test_Reaction_Yes);
        
    end
        
end
FDR_PHB = mafdr(pvalue_PHB,'BHFDR','true');


% MDF
[All_Directions,Nr_EFMs] = size(EFM_Count_Split_MDF);
used = 0;
for l=1:All_Directions  % Loop through all reactions (including reverse)
    clear Test_Reaction
    clear Test_Reaction_Null
    clear Test_Reaction_Yes
    if abs(sum(EFM_Count_Split_MDF(l,:))) ~= Nr_EFMs && ~sum(EFM_Count_Split_MDF(l,:)) == 0 % if the reaction is not used in all EFMs and is used (e.g. it's not zero)
        used = used + 1;
        Test_Reaction = EFM_Count_Split_MDF(l,:); % Take the reaction to be tested
        Reac_MDF{used,1} = Reaction_Name_Before_MDF(l);
        Score_MDF(used,1) = MDF_Score_trans(l);
        for q = 1:length(Test_Reaction) 
            if  Test_Reaction(q) == 1
                Test_Reaction_Null(q) = 0;
                Test_Reaction_Yes(q) = MDF_total(q);
            else
                Test_Reaction_Null(q) = MDF_total(q);
                Test_Reaction_Yes(q) = 0;
            end
        end
        Test_Reaction_Null(Test_Reaction_Null==0) = [];
        Test_Reaction_Yes(Test_Reaction_Yes==0) = [];
        [pvalue_MDF(used,1),h_MDF(used)] = ranksum(Test_Reaction_Null,Test_Reaction_Yes);
        
    end
        
end
FDR_MDF = mafdr(pvalue_MDF,'BHFDR','true');

