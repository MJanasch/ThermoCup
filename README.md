![alt text](ThermoCup_Logo.png "Thermodynamics in Cupriavidus necator aka ThermoCup")

# ThermoCup - Thermodynamic analysis of the core metabolism in Cupriavidus necator in different substrate utilization scenarios

Code for computational analysis used in Janasch et al., 2021

---

### 1. EFM Enumeration, filtering and thermodynamic analysis

##### Scripts

'ThermoCup_EnumerateEFMs.m' - script to enumerate EFMs (Elementary flux modes) and perform MDF analysis, but no variability analysis
The substrate needs to be defined ("Frc", "Suc" or "For"), as is described in the script.
Requires Support scripts:
'Support/ThermoCup_ReadExcel.m'
'Support/ThermoCup_MDF.m'
'Support/ThermoCup_SolveLP_MDF.m'
'Support/ThermoCup_CalculateTransportEnergy.m'
'Support/ThermoCup_FindIndex.m'


##### Input data: 
'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure


##### Output:

Yield values for EFMs leading to either PHB and/or biomass production from the corresponding substrate, as well as the corresponding MDF values.
The values for each substrate were saved as 'Data/ThermoCup_InputDataforScores_Frc.mat', 'ThermoCup_InputDataforScores_Suc.mat' and 'ThermoCup_InputDataforScores_For.mat', to be used for "2. Reaction Usage" and "3. Yield and MDF Scores"

Recommended:
Save all enumerated EFMs for all three substrates ("mnet_Frc", "mnet_Suc", "mnet_For"), in one file containing all enumerated EFMs 'ThermoCup_All_Enumerated_EFMs.mat' into the Data subfolder, to skip the EFMTool enumeration in future applications.

---

### 2. Reaction Usage

##### Scripts:
'ThermoCup_GetReactionUsage.m' - script to analyze reaction usage, needs to be run independently for each strain and substrate. Splits reversible reactions and calculates number of EFMs that use a particular reaction, as a fraction of all EFMs analyzed. Data is used for Supplementary figure S2)


##### Input data:
'Data/ThermoCup_All_Enumerated_EFMs.mat' - File containing all enumerated EFMs, for all three substrates ("mnet_Frc", "mnet_Suc", "mnet_For"), created in Step 1.

EFM_Count data as part of 'Data/ThermoCup_InputDataforScores_Frc.mat', 'ThermoCup_InputDataforScores_Suc.mat' and 'ThermoCup_InputDataforScores_For.mat' for the three substrates, respectively.

##### Output:

 `EFM_Count_Sum` as a vector containing the percentage reaction usage for each reaction.

---

### 3. Production envelope calculation via Flux Balance Analysis (FBA)

##### Scripts:
Scripts for FBA are found in the folder 'CupriavidusCoreFBA'.
'ThermoCup_FBA.m' - script to calculate maximum yields as well as production envelope. Reaction constraints need to be adjusted according to substrate utilization.

Requires 'ThermoCup_ReadExcel_FBA.m' - function to read in the CupriaCore FBA model.
Requires furthermore installed COBRA toolbox, see https://opencobra.github.io/cobratoolbox/stable/.


##### Input data:
'Model/Cupriavidus_Core_ALL_FBA.xlsx' - excel file containing the FBA model structure

##### Output:
Three "Result" tables (one for each strain) containing data points for Biomass and PHB yield, forming the production envelope. Needs to be run for each substrate.


---

### 4. Yield and MDF Scores

##### Scripts:
'ThermoCup_GetReactionScores.m' - script to calculate biomass and PHB yield scores and MDF scores, dependent on the substrate. Also calculates significance via Wilcoxon Rank Sum test and FDR.

##### Input data:
'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure

'Data/ThermoCup_All_Enumerated_EFMs.mat' - File containing all enumerated EFMs, for all three substrates ("mnet_Frc", "mnet_Suc", "mnet_For"), created in Step 1.

Biomass/PHB yield data and MDF data as part of 'Data/ThermoCup_InputDataforScores_Frc.mat', 'ThermoCup_InputDataforScores_Suc.mat' and 'ThermoCup_InputDataforScores_For.mat' for the three substrates, respectively.

##### Output:

 Yield score data for Biomass and PHB, MDF score data, corresponding P-values and FDR, to be used in Figure 3, shown in Supplementary table TS2-4.



---

### 5. Gibbs free energy change variability analysis for the subnetwork in Figure 3

##### Scripts:
'ThermoCup_SmallNetworkMDF.m' - script to calculate Gibbs free emergy change variability of reactions of the canonical TCA and the reverse MDH alternative.

Requires Support scripts:
'Support/ThermoCup_SolveLP_ERA_MDF_SmallNetwork.m'
'Support/ThermoCup_MDF_SmallNetwork.m'
'Support/ThermoCup_MDF_EVA_dG.m'
'Support/ThermoCup_FindIndex.m'



##### Input data:
'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure 

'Data/ThermoCup_All_Enumerated_EFMs.mat' - File containing all enumerated EFMs, for all three substrates ("mnet_Frc", "mnet_Suc", "mnet_For"), here "mnet_For" was used.


##### Output:

Two datastructures ("DF_canonTCA" and "DF_revMDH") containing the min/max Gibbs free energy change for the reactions involved in either exemplary case.

---

### 6. Metabolite variability analysis and metabolite sampling preparation

##### Scripts:
'ThermoCup_MDF_EVA_DataOut.m' - script similar to 'ThermoCup_EnumerateEFMs.m', with addition of the variability analysis. Loops through all EFMs, determines MDF and performs metabolite variability analysis and saves all data in .txt files to be used as input for metabolite sampling and for 'ThermoCup_AS2_TransformMetRanges.m'. Needs to be run for each substrate seperately.

Requires Support scripts:
'Support/ThermoCup_ReadExcel.m'
'Support/ThermoCup_MDF.m'
'Support/ThermoCup_SolveLP_MDF.m'
'Support/ThermoCup_CalculateTransportEnergy.m'
'Support/ThermoCup_FindIndex.m'
'Support/ThermoCup_MDF_EVA.m'

##### Input:
'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure

'Data/ThermoCup_All_Enumerated_EFMs.mat' - File containing all enumerated EFMs, for all three substrates ("mnet_Frc", "mnet_Suc", "mnet_For"), created in Step 1.


##### Output:
Depending on the substrate chosen, a folder will be created ("Output_for_Sampling_Frc", "Output_for_Sampling_Suc", "Output_for_Sampling_For"). Each folder contains subfolder, containing a textfile for each EFM containing information to be used as input for metabolite sampling and subsequent analysis:
- Initial Concentration ("Conc_Init")
- standard Gibbs Free Energy change ("dG0")
- MDF ("MDF")
- Metabolite Ranges from variability analysis ("Met_Ranges")
- Names of reactions and metabolites included in the specific EFM ("Name_References")
- Matrix for cofactor ratios ("Ratio_Matrix")
- Stoichiometric Matrix ("S_Matrix")


---

### 7. Random sampling of metabolite concentration/driving force

##### Scripts:
'Metabolite_sampling/ThermoCup_HnR_MetaboliteSampling.py' - python script that reads in data from the folders created under  5. for each EFM, and performs hit-and-run metabolite sampling, 5000 steps. Needs to be run for all EFMs for each substrate separately.

##### Input:
Data created for each EFM under 5.
- Initial Concentration ("Conc_Init")
- standard Gibbs Free Energy change ("dG0")
- MDF ("MDF")
- Metabolite Ranges from variability analysis ("Met_Ranges")
- Names of reactions and metabolites included in the specific EFM ("Name_References")
- Matrix for cofactor ratios ("Ratio_Matrix")
- Stoichiometric Matrix ("S_Matrix")


##### Output:
For each EFM, two text-files are being created. They contain the median driving force for each reaction in the corresponding EFM (from 5000 steps) and their MAD, respectively. The text-files are saved in the folders "Substrat_PHB_Median" and "Substrate_PHB_MAD".

Optional:
Output of all 5000 metabolite concentration sets for each EFM.

---

### 8. Analysis of sampling results

##### Scripts:
AS - "After Sampling"

1) 'ThermoCup_AS1_SortnPrepareDF' - script that loads driving forces, MDF and reaction names for each EFM and puts into more usable data structure in matlab. Furthermore it sorts the EFMs after 1) their MDF value and 2) EFMs with the same MDF value are sorted by the mean of their driving forces. Needs to be run sperataly for each substrate.


2) 'ThermoCup_AS2_TransformMetRanges.m' - script that reads in metabolite ranges from the variability analysis and sorted indices from 1), to sort metabolite ranges and put data into structure and outputs text-files to create supplementary figures S4-S6. Needs to be run sperataly for each substrate.


3) 'ThermoCup_AS3_GetDFNorm_and_DFSumsAcCoA.m' - script that reads in the sorted DF-datastructure from 1), relates each reaction's DF to the MDF of the EFM, thereby finding how far away the DF is from the MDF. This identifies which reactions are at the thermodynamic limit of the MDF. In a second part the script sums all DFs leading to acetyl-CoA. The script creates the data used for figure 4, as well as supplementary figures S7-S10. Needs to be run sperataly for each substrate.

4) 'ThermoCup_AS4_GetDF_Map.m' - script to create DF overview for the two metabolic strategies found for formate.


##### Input:
1) From random sampling and variability analysis:
- Name_References (Names of all reactions and metabolites for each EFM)
- MDF (MDF value for each EFM)
- PHB_Median (Median DFs for each reaction for each EFM)
- PHB_MAD (MADs of the DFs for each reaction for each EFM)


2) 
- 'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure (Supplementary file S1 in the article)
From random sampling and variability analysis:
- Name_References (Names of all reactions and metabolites for each EFM)
- MDF (MDF value for each EFM)
- Met_Ranges (Metabolite Ranges from vriability analysis)

- Sorted indices from 1) ('Indices_DoubleSorted_Frc.mat', 'Indices_DoubleSorted_Suc.mat' and 'Indices_DoubleSorted_For.mat')


3)
- Sorted DF from 1) ('Sorted_EFMs_Frc_DF.mat', 'Sorted_EFMs_Suc_DF.mat', 'Sorted_EFMs_For_DF.mat')


4) 
- 'Model/Cupriavidus_Core_XII_PHB.xlsx' - excel file containing the model structure (Supplementary file S1 in the article)
- Sorted DF from 1) for formate ('Sorted_EFMs_For_DF.mat')

##### Output:
1) Datastructure containing DFs for each reaction of each EFM, EFMs sorted by 1) MDF and 2) mean over all DFs. Saved as 'Sorted_EFMs_Frc_DF.mat', 'Sorted_EFMs_Suc_DF.mat', 'Sorted_EFMs_For_DF.mat'.
Additionally saves indices of sorted EFMs to be able to sort metabolite ranges accordingly, as 'Indices_DoubleSorted_Frc.mat', 'Indices_DoubleSorted_Suc.mat' and 'Indices_DoubleSorted_For.mat'.


2) Textfiles containing Max/min concentration for each metabolite for each EFM, sorted by 1) MDF and 2) mean DF. If a metabolite is not part of the corresponding EFM, the value was arbitrarily set to 1000. The header needed for the R script is created as well.

3) Matrices containing DF and normalized DF for each reaction for each EFM in "DF_for_Plot_Matrix" and "DF_for_Plot_Matrix_Norm", respectively.

In the second part the sum of all DFs towards acetyl-CoA for each strains is found as "Sources_wt", "Sources_xfpk" and "Sources_aclAB" (no aclAB for formate).

4) Outputs two datastructures "Stat_DF_1" and "Stat_DF_2", containing mean DF for each reaction over all EFMs participating in the corresponding metabolic strategy.


---


## Note

- R Scripts to create the figures used in the manuscript can be provided by the author upon request.


---

## Dependencies

EFMtool for Matlab (https://csb.ethz.ch/tools/software/efmtool.html)

Mathworks Matlab (Performed on MATLAB_R2018a)

Gurobi solver for FBA, free academic license (https://www.gurobi.com)

Python 3 (3.8.5) with the following libraries:
- sys, os, numpy, time, math, scipy


---

## Authors

Markus Janasch, SciLifeLab/KTH (markus.janasch@scilifelab.se)