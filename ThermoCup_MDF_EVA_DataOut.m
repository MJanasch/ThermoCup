%% Combine MDF, EVA and Data Output for Sampling
% Perform EFM filtering, MDF calculation as well as EVA, save output in
% text-files to be read into Hit-n-Run metabolite sampling


%% Suppress warnings
warning('off','all') % Warnings will appear because the parent folder "Output_for_Sampling" exists after it's created on the first occasion

%% Load model
addpath('Model/')
modelFile = '/Model/Cupriavidus_Core_XII_PHB.xlsx';
CupriaCore = ThermoCup_ReadExcel(modelFile);

%% Load enumerated EFMs
addpath('Data/')
load('ThermoCup_All_Enumerated_EFMs.mat');

% Set Folder
currentFolder = pwd;

%% Check for biomass and PHB production

Substrate = 'Frc';
% Substrate = 'Suc';
% Substrate = 'For';

if exist('Substrate','var') == 1
    if strcmp(Substrate,'Frc')
        mnet = mnet_Frc;
        MolWeight = 0.18;     % Molecular weight of fructose
        Xfpk1 = 72;
        Xfpk2 = 73;
        AclAB = 74;
        Model_Variant = CupriaCore.Variant_Frc;
        disp('Substrate is Frc')
    elseif strcmp(Substrate,'Suc')
        mnet = mnet_Suc;
        MolWeight = 0.1181;   % Molecular weight of succinate
        Xfpk1 = 72;
        Xfpk2 = 73;
        AclAB = 74;
        Model_Variant = CupriaCore.Variant_Suc;
        disp('Substrate is Suc')
    elseif strcmp(Substrate,'For')
        mnet = mnet_For;
        MolWeight = 0.046;      % Molecular weight of formate
        Xfpk1 = 73;
        Xfpk2 = 74;
        AclAB = 75;
        Model_Variant = CupriaCore.Variant_For;
        disp('Substrate is For')
    else
        disp('Choose a right substrate')
    end
else
    disp('Define a substrate, please choose between "Frc", "Suc" and "For".')
end
Nr_Rxns = length(mnet.reactionNames);


%% Select EFMs for reaction usage
% mnet.efms(4,i) is always biomass
% mnet.efms(9,i) is always PHB exchange

tic
r=0;
r_yes = 0;
for i = 1:length(mnet.efms)
    if mnet.efms(9,i) ~= 0 && ((mnet.efms(74,i) == 0 && mnet.efms(72,i) == 0 && mnet.efms(73,i) == 0) || ... % WT
            (mnet.efms(74,i) ~= 0 && mnet.efms(72,i) == 0 && mnet.efms(73,i) == 0) || ...   % aclAB
            (mnet.efms(74,i) == 0 && ...
            ((mnet.efms(72,i) ~= 0 && mnet.efms(73,i) ~= 0) ||...
            (mnet.efms(72,i) ~= 0 && mnet.efms(73,i) == 0) || ...
            (mnet.efms(72,i) == 0 && mnet.efms(73,i) ~= 0))))
        r=r+1;
        Result(i) = 1;
        Index(r) = i;
        [MDF(r),modelOut] = ThermoCup_MDF(Model_Variant,CupriaCore.metabolites,CupriaCore.cofaratios,mnet.efms(:,i));
        if MDF(r) > 0    
            r_yes = r_yes +1;
            EVA(r_yes).Results = ThermoCup_MDF_EVA(modelOut,MDF(r),0.9);
            %EVA_dG(r_yes).Results = ThermoCup_MDF_EVA_dG(modelOut,MDF(r),0.9);
            MDF_yes = r_yes;
            model_temp = modelOut;
        
                %% S-matrix

                S_temp = model_temp.S_MOD;
                Strnsp_temp = model_temp.Strnsp_MOD;
                [n,m] = size(S_temp);
               
                
                
%                 RunVar = m:-1:1;
%                 for b = 1:length(RunVar)
%                     if strncmp(model_temp.reactions.rxns(RunVar(b)),'Ex_',3)
%                         S_temp(:,RunVar(b))            = [];
%                         Strnsp_temp(:,RunVar(b))       = [];
%                         
%                     elseif strcmp(model_temp.reactions.rxns(RunVar(b)),'BIOMASS')
%                         S_temp(:,RunVar(b))            = [];
%                         Strnsp_temp(:,RunVar(b))       = [];
%                     elseif strncmp(model_temp.reactions.rxns(RunVar(b)),'Trp_ForO',7)
%                         S_temp(:,RunVar(b))            = [];
%                         Strnsp_temp(:,RunVar(b))       = [];
%                     elseif strcmp(model_temp.reactions.rxns(RunVar(b)),'PHB synthase')
%                         S_temp(:,RunVar(b))            = [];
%                         Strnsp_temp(:,RunVar(b))       = [];
%                     end
%                 end

                % Save S-matrix into textfile
                mkdir Output_for_Sampling S_Matrix
                cd('Output_for_Sampling/S_Matrix/');
                filename = ['S_Matrix_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');


                for b = 1:size(S_temp,1)
                    fprintf(fileID,'[');
                    fprintf(fileID,' %g',S_temp(b,1:end-1));
                    fprintf(fileID,' %g],\n',S_temp(b,end));
                end
                %fprintf(fileID,']');
                fclose(fileID);
                cd(currentFolder)


                %% dG0's
                [n_dG,m_dG] = size(S_temp);
                %for Rxn_Index = 1:m_dG
                dG0 = transpose(S_temp)*model_temp.metabolites.dfG + transpose(Strnsp_temp)*model_temp.reactions.dtG;
                %end

                % Save dG0 into textfile
                dG0=transpose(dG0);

                mkdir Output_for_Sampling dG0
                cd('Output_for_Sampling/dG0/');
                filename = ['dG0_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');
                fprintf(fileID,'[');
                fprintf(fileID,' %g,',dG0(1,1:end-1));
                fprintf(fileID,' %g]',dG0(1,end));
                %fprintf(fileID,']');
                fclose(fileID);
                cd(currentFolder)

                %% Metabolite Ranges
                EVA_temp = EVA(MDF_yes).Results;
                numtot_Mets_temp = length(EVA_temp);
                for Index_Met = 1:numtot_Mets_temp
                    Conc_Ranges(Index_Met,1) = EVA_temp(Index_Met,1).Conc_Out;
                    Conc_Ranges(Index_Met,2) = EVA_temp(Index_Met,2).Conc_Out;
                end


                %Save metabolite ranges into textfile
                mkdir Output_for_Sampling Met_Ranges
                cd('Output_for_Sampling/Met_Ranges/');
                filename = ['MetRange_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');
                %fprintf(fileID,'[');
                for Index_Met = 1:numtot_Mets_temp
                    fprintf(fileID,'[%g %g],\n',Conc_Ranges(Index_Met,1),Conc_Ranges(Index_Met,2));
                end
                %fprintf(fileID,']');
                fclose(fileID);
                cd(currentFolder)

                %% Ratio Matrix
                % 4 Ratios
                Nr_Ratios = 4;
                Ratio_Matrix = zeros(numtot_Mets_temp,Nr_Ratios);
                for j = 1:Nr_Ratios
                    Ratio_Split = split(model_temp.cofaratios.ratios(j),"/");
                    Over = Ratio_Split(1);
                    Under = Ratio_Split(2);
                    Ratio_Matrix(MJanasch_FindIndex(model_temp.metabolites.mets,Over),j) = 1;
                    Ratio_Matrix(MJanasch_FindIndex(model_temp.metabolites.mets,Under),j) = -1;
                end


                % Save R-matrix into textfile
                mkdir Output_for_Sampling Ratio_Matrix
                cd('Output_for_Sampling/Ratio_Matrix/');
                filename = ['Ratio_Matrix_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');


                for b = 1:size(Ratio_Matrix,1)
                    fprintf(fileID,'[');
                    fprintf(fileID,' %g',Ratio_Matrix(b,1:end-1));
                    fprintf(fileID,' %g],\n',Ratio_Matrix(b,end));
                end
                fclose(fileID);
                cd(currentFolder)

                %% Actual MDF value for each EFM

                % Save MDF into textfile
                mkdir Output_for_Sampling MDF
                cd('Output_for_Sampling/MDF/');
                filename = ['MDF_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');
                fprintf(fileID,'%g',MDF(r));
                fclose(fileID);
                cd(currentFolder)



                %% Save initial concentration set

                conc_temp = model_temp.metabolites.conc;
                
                
                
                
                % Save MDF into textfile
                mkdir Output_for_Sampling Conc_Init
                cd('Output_for_Sampling/Conc_Init/');
                filename = ['Conc_Init_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');
                fprintf(fileID,'[');
                fprintf(fileID,' %g,',conc_temp(1:end-1,1));
                fprintf(fileID,' %g]',conc_temp(end,1));
                fclose(fileID);
                cd(currentFolder)
            
                %% Save names of reactions and metabolites involved
                Rxn_Names = model_temp.reactions.rxns_dG;
                
                for p=1:m
                    Rxn_Names{p} = ['R_',Rxn_Names{p}];
                end
                
                Met_Names = model_temp.metabolites.mets;
                for q=1:n
                    Met_Names{q} = ['M_',Met_Names{q}];
                end
                
                mkdir Output_for_Sampling Name_References
                cd('Output_for_Sampling/Name_References/');
                filename = ['Name_References_EFM_Nr_',int2str(MDF_yes),'_',Substrate,'.txt'];
                fileID = fopen(filename,'w');
                Name_File_Content = [Rxn_Names;Met_Names];
                NFC = Name_File_Content.';
                fprintf(fileID,'%s\n',Name_File_Content{:});
                fclose(fileID);
                cd(currentFolder)
        end
        yield(r) = mnet.efms(4,i)/(mnet.efms(1,i)*MolWeight);
        for o = 1:Nr_Rxns
            if mnet.efms(o,i) < 0
                EFM_Count(o,r) = -1;
            elseif mnet.efms(o,i) > 0
                EFM_Count(o,r) = 1;
            else
                EFM_Count(o,r) = 0;
            end
        end
        if rem(r,1000) == 0
            disp(['Number of EFMs found: ',num2str(r)])
            toc
        end
    else
        Result(i) = 0;
    end
    
end
disp(['Total number of EFMs found: ',num2str(r)])
toc
