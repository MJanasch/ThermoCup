%% Calculate Change in Gibbs Free Energy for transport processes
% Markus Janasch
% Created: 20-08-27     Last Modified: -



function [dtG_Names, dtG_Values] = ThermoCup_CalculateTransportEnergy(RT,pH_in,pH_out,membrane_potential,F)

conc_proton_in    = 10^-pH_in;
conc_proton_out   = 10^-pH_out;

%membrane_potential = -0.15; % [V]

%F = 96.48533212; % [C/mmol]

% Fructose
dfG_Fructose_0_in = -915.5 - 12*(RT*log(conc_proton_in));
dfG_Fructose_1_in = -857.7 - 11*(RT*log(conc_proton_in));
dfG_Fructose_2_in = -785.1 - 10*(RT*log(conc_proton_in));
dfG_Fructose_3_in = -705.8 - 9*(RT*log(conc_proton_in));

dfG_Fructose_0_out = -915.5 - 12*(RT*log(conc_proton_out));
dfG_Fructose_1_out = -857.7 - 11*(RT*log(conc_proton_out));
dfG_Fructose_2_out = -785.1 - 10*(RT*log(conc_proton_out));
dfG_Fructose_3_out = -705.8 - 9*(RT*log(conc_proton_out));

dfG_Fructose_in = -RT*log(exp(-dfG_Fructose_0_in/RT)+exp(-dfG_Fructose_1_in/RT)+exp(-dfG_Fructose_2_in/RT)+exp(-dfG_Fructose_3_in/RT));
dfG_Fructose_out= -RT*log(exp(-dfG_Fructose_0_out/RT)+exp(-dfG_Fructose_1_out/RT)+exp(-dfG_Fructose_2_out/RT)+exp(-dfG_Fructose_3_out/RT));


% CO2
dfG_CO2_0_in = -623.1 - 2*(RT*log(conc_proton_in));
dfG_CO2_1_in = -586.8 - 1*(RT*log(conc_proton_in));
dfG_CO2_2_in = -527.8 - 0*(RT*log(conc_proton_in));

dfG_CO2_0_out = -623.1 - 2*(RT*log(conc_proton_out));
dfG_CO2_1_out = -586.8 - 1*(RT*log(conc_proton_out));
dfG_CO2_2_out = -527.8 - 0*(RT*log(conc_proton_out));

dfG_CO2_in = -RT*log(exp(-dfG_CO2_0_in/RT)+exp(-dfG_CO2_1_in/RT)+exp(-dfG_CO2_2_in/RT));
dfG_CO2_out= -RT*log(exp(-dfG_CO2_0_out/RT)+exp(-dfG_CO2_1_out/RT)+exp(-dfG_CO2_2_out/RT));


% EtOH
dfG_EtOH_0_in = -181.6 - 6*(RT*log(conc_proton_in));
dfG_EtOH_0_out = -181.6 - 6*(RT*log(conc_proton_out));

dfG_EtOH_in = -RT*log(exp(-dfG_EtOH_0_in/RT));
dfG_EtOH_out= -RT*log(exp(-dfG_EtOH_0_out/RT));


% Pyruvate
dfG_Pyruvate_0_in = -485.8 -4*(RT*log(conc_proton_in));
dfG_Pyruvate_1_in = -472.3 -3*(RT*log(conc_proton_in));

dfG_Pyruvate_0_out = -485.8 -4*(RT*log(conc_proton_out));
dfG_Pyruvate_1_out = -472.3 -3*(RT*log(conc_proton_out));

dfG_Pyruvate_in = -RT*log(exp(-dfG_Pyruvate_0_in/RT)+exp(-dfG_Pyruvate_1_in/RT));
dfG_Pyruvate_out= -RT*log(exp(-dfG_Pyruvate_0_out/RT)+exp(-dfG_Pyruvate_1_out/RT));


% Phosphate
dfG_Phosphate_0_in = -1143.5 - 3*(RT*log(conc_proton_in));
dfG_Phosphate_1_in = -1137.3 - 2*(RT*log(conc_proton_in));
dfG_Phosphate_2_in = -1096.1 - 1*(RT*log(conc_proton_in));
dfG_Phosphate_3_in = -1020.0 - 0*(RT*log(conc_proton_in));

dfG_Phosphate_0_out = -1143.5 - 3*(RT*log(conc_proton_out));
dfG_Phosphate_1_out = -1137.3 - 2*(RT*log(conc_proton_out));
dfG_Phosphate_2_out = -1096.1 - 1*(RT*log(conc_proton_out));
dfG_Phosphate_3_out = -1020.0 - 0*(RT*log(conc_proton_out));

dfG_Phosphate_in = -RT*log(exp(-dfG_Phosphate_0_in/RT)+exp(-dfG_Phosphate_1_in/RT)+exp(-dfG_Phosphate_2_in/RT)+exp(-dfG_Phosphate_3_in/RT));
dfG_Phosphate_out= -RT*log(exp(-dfG_Phosphate_0_out/RT)+exp(-dfG_Phosphate_1_out/RT)+exp(-dfG_Phosphate_2_out/RT)+exp(-dfG_Phosphate_3_out/RT));


% Formate
dfG_Formate_0_in = -372.1 -2*(RT*log(conc_proton_in));
dfG_Formate_1_in = -351.0 -1*(RT*log(conc_proton_in));

dfG_Formate_0_out = -372.1 -2*(RT*log(conc_proton_out));
dfG_Formate_1_out = -351.0 -1*(RT*log(conc_proton_out));

dfG_Formate_in = -RT*log(exp(-dfG_Formate_0_in/RT)+exp(-dfG_Formate_1_in/RT));
dfG_Formate_out= -RT*log(exp(-dfG_Formate_0_out/RT)+exp(-dfG_Formate_1_out/RT));


% Succinate
dfG_Succinate_0_in = -746.6 - 6*(RT*log(conc_proton_in));
dfG_Succinate_1_in = -722.6 - 5*(RT*log(conc_proton_in));
dfG_Succinate_2_in = -690.4 - 4*(RT*log(conc_proton_in));

dfG_Succinate_0_out = -746.6 - 6*(RT*log(conc_proton_out));
dfG_Succinate_1_out = -722.6 - 5*(RT*log(conc_proton_out));
dfG_Succinate_2_out = -690.4 - 4*(RT*log(conc_proton_out));

dfG_Succinate_in = -RT*log(exp(-dfG_Succinate_0_in/RT)+exp(-dfG_Succinate_1_in/RT)+exp(-dfG_Succinate_2_in/RT));
dfG_Succinate_out= -RT*log(exp(-dfG_Succinate_0_out/RT)+exp(-dfG_Succinate_1_out/RT)+exp(-dfG_Succinate_2_out/RT));


% NH3 Ammonia
dfG_NH3_0_in = -26.5 -3*(RT*log(conc_proton_in));
dfG_NH3_pos1_in = -79.3 -4*(RT*log(conc_proton_in));

dfG_NH3_0_out = -26.5 -3*(RT*log(conc_proton_out));
dfG_NH3_pos1_out = -79.3 -4*(RT*log(conc_proton_out));

dfG_NH3_in = -RT*log(exp(-dfG_NH3_0_in/RT)+exp(-dfG_NH3_pos1_in/RT));
dfG_NH3_out= -RT*log(exp(-dfG_NH3_0_out/RT)+exp(-dfG_NH3_pos1_out/RT));


% O2
dfG_O2_in  = 16.4;
dfG_O2_out = 16.4;


% H2O
dfG_H2O_in  = -237.2 -2*(RT*log(conc_proton_in));
dfG_H2O_out = -237.2 -2*(RT*log(conc_proton_out));


% Proton
dfG0_H_in  = - 1*(0+RT*log(conc_proton_in));
dfG0_H_out = - 1*(0+RT*log(conc_proton_out));

%% dG associated with transport (transport of charge, protons, H-atoms)

% for Fructose-abc: no proton transport, no charge, 12 protons in Frc
dtG_Fructose = -12*RT*log(conc_proton_out)+12*RT*log(conc_proton_in); % [J/mmol]


% for Succinate-abc: 2 negative charges, no protons, 4 hydrogen atoms on succinate2-
dtG_Succinate_ABC = -4*RT*log(conc_proton_out)+4*RT*log(conc_proton_in) + F*membrane_potential*-2; % [J/mmol]


% for Succinate proton symport: no charge, 2 protons, 4 hydrogen atoms on succinate2-
dtG_Succinate_SYM = 2*dfG0_H_in - 2*dfG0_H_out -6*RT*log(conc_proton_out)+6*RT*log(conc_proton_in); % [J/mmol]


% Formate dehydrogenase, membrane bound: no protons, no hydrogen atoms, not even molecules, just 2 electrons (i.e. 2 negative charges)!
dtG_Formate_MBRN = F*membrane_potential*-2; % [J/mmol]


% Formate symport: 1 proton symport, no charge, 1 hydrogen atom on formate
dtG_Formate_SYM = 1*dfG0_H_in - 1*dfG0_H_out -2*RT*log(conc_proton_out)+2*RT*log(conc_proton_in); % [J/mmol]


% for CO2: proton symport, no charge because leveled out, 1 proton transported by itself, one proton on CO2 (as HCO3^-1)
dtG_CO2 = -1*dfG0_H_in + 1*dfG0_H_out + -2*RT*log(conc_proton_in) + 2*RT*log(conc_proton_out); % [J/mmol], stoichiometry of 2 for the last two summands because 1 proton transported attached to CO2tot, one as the symport, also quite close to what Niebel et al had


% for H2O: free, unhindered "diffusion"
dtG_H2O = 0;


% for O2: free, unhindered "diffusion"
dtG_O2 = 0;


% for Phosphate-abc: no proton transport, 2 negative charges, 1 proton on phosphate
dtG_Phosphate_ABC = -1*RT*log(conc_proton_out)+1*RT*log(conc_proton_in) + F*membrane_potential*-2; % [J/mmol], 2 negative charges transported on Pi, one hydrogen on Pi2-


% for Phosphate proton symport: 2 proton transported, no charges, 1 proton on phosphate
dtG_Phosphate_SYM = -2*dfG0_H_out + 2*dfG0_H_in + -3*RT*log(conc_proton_out)+3*RT*log(conc_proton_in); % [J/mmol], 2 negative charges transported on Pi, one hydrogen on Pi2-


% for Pyruvate: proton symport, no charge because leveled out, 1 proton transported by itself, 3 protons on pyruvate (as C3H3O3^1-)
dtG_Pyruvate = -1*dfG0_H_in + 1*dfG0_H_out + -4*RT*log(conc_proton_in) + 4*RT*log(conc_proton_out); % [J/mmol]


% for NH3 Ammonia: proton ANTIPORT, no charge, 1 proton transported OUT of the cell, 4 hydrogens transported inside via NH4+, but one proton going out, therefore only stoichiometry 3
dtG_NH3 = -1*dfG0_H_in + 1*dfG0_H_out + -3*RT*log(conc_proton_out) + 3*RT*log(conc_proton_in); % [J/mmol]



% for proton-translocating transhydrogenase: 2 protons transported out, 2 positive charges, (NADPH[c] + NAD[c] + 2 H+[c] <=> NADP[c] + NADH[c] + 2 H+[e])
%dtG_THD = 2*dfG0_H_in + -2*dfG0_H_out + 2*RT*log(conc_proton_out) + -2*RT*log(conc_proton_in) + F*membrane_potential*2; % [J/mmol]



%% Summarizing all transport dG's
dtG_Names = {'Trp_FrcABC';
    'Trp_O2';
    'Trp_CO2';
    'Trp_H2O';
    'Trp_piABC';
    'Trp_pi';
    'Trp_Pyr';
    'Trp_NH3';
    'Trp_SucABC';
    'Trp_SucSYM';
    'Trp_FDH_Mem';
    'Trp_ForSYM'};

dtG_Values = [dtG_Fructose;
    dtG_O2;
    dtG_CO2;
    dtG_H2O;
    dtG_Phosphate_ABC;
    dtG_Phosphate_SYM;
    dtG_Pyruvate;
    dtG_NH3;
    dtG_Succinate_ABC;
    dtG_Succinate_SYM;
    dtG_Formate_MBRN;
    dtG_Formate_SYM];

end