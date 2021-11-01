## Hit-and-Run Sampling, adapted from Johannes Asplund-Samuelsson (https://github.com/Asplund-Samuelsson)

# Import libraries
import sys, os
import numpy as np
import time
import math
from scipy import stats


#######################################################################################################
## Names of input files and output files need to be changed according to which substrate is being used!
#######################################################################################################




EFM_Nr = sys.argv[1]


#########---Read in Data---#########

###-----Load Stoichiometric Matrix-----###

S_Matrix_file_name = sys.argv[2]
#S_Matrix_file_name = "/S_Matrix/S_Matrix_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+S_Matrix_file_name
S_Matrix_file = open(path,"r+")
S_Matrix_file_contents = S_Matrix_file.read()

S_Matrix_file_contents = S_Matrix_file_contents[:-2]
S_Matrix_file_contents = S_Matrix_file_contents.replace("\n"," ")
S_Matrix_file_contents = S_Matrix_file_contents.split(", ")

S_Matrix = []
for line in S_Matrix_file_contents:
    line = line[1:-1]
    line = list(line.split(" "))
    line = line[1:]
    
    line_float = [float(entry) for entry in line]
    S_Matrix.append(line_float)

S_Matrix = np.array(S_Matrix)
#print(S_Matrix)


###-----Load Standard Change of Gibbs Free Energy Values-----###
dG0_file_name = sys.argv[3]
#dG0_file_name = "/dG0/dG0_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+dG0_file_name
dG0_file = open(path,"r+")
dG0_file_contents = dG0_file.read()
dG0_file_contents = dG0_file_contents[2:-1]
dG0_file_contents = dG0_file_contents.split(', ')
dG0_float = [float(entry) for entry in dG0_file_contents]

dG0 = np.array(dG0_float)

# RT is a constant
T=303.15
R=8.3145e-3
RT = R*T



###-----Load Metabolite Concentration ranges-----###
MetRange_file_name = sys.argv[4]
#MetRange_file_name = "/Met_Ranges/MetRange_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+MetRange_file_name
MetRange_file = open(path,"r+")
MetRange_file_contents = MetRange_file.read()
MetRange_file_contents = MetRange_file_contents[:-2]
MetRange_file_contents = MetRange_file_contents.replace("\n"," ")
MetRange_file_contents = MetRange_file_contents.split(", ")

MetRange = []
for line in MetRange_file_contents:
    line = line[1:-1]
    
    line = list(line.split(" "))
    
    line_float = [float(entry)/1000 for entry in line]
    MetRange.append(line_float)

#MetRange = np.log(np.array(MetRange))
MetRange = np.round(np.log(np.array(MetRange)),3)

#print(MetRange)


###-----Load MDF Value-----###
MDF_file_name = sys.argv[5]
#MDF_file_name = "/MDF/MDF_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+MDF_file_name
MDF_file = open(path,"r+")
MDF_file_contents = MDF_file.read()
#MDF = round(float(MDF_file_contents),2)
MDF = float(MDF_file_contents)

###-----Load Starting Concentration set-----###
Conc_Init_file_name = sys.argv[6]
#Conc_Init_file_name = "/Conc_Init/Conc_Init_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+Conc_Init_file_name
Conc_Init_file = open(path,"r+")
Conc_Init_file_contents = Conc_Init_file.read()
Conc_Init_file_contents = Conc_Init_file_contents[2:-1]
Conc_Init_file_contents = Conc_Init_file_contents.split(', ')
Conc_Init_float = [float(entry) for entry in Conc_Init_file_contents]
c_0 = np.round(np.log(np.array(Conc_Init_float)),3)

#c_0 = np.log(np.array(Conc_Init_float))
#print(c_0)



###-----Load Ratio Matrix-----###
R_Matrix_file_name = sys.argv[7]
#R_Matrix_file_name = "/Ratio_Matrix/Ratio_Matrix_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+R_Matrix_file_name
R_Matrix_file = open(path,"r+")
R_Matrix_file_contents = R_Matrix_file.read()

R_Matrix_file_contents = R_Matrix_file_contents[:-2]
R_Matrix_file_contents = R_Matrix_file_contents.replace("\n"," ")
R_Matrix_file_contents = R_Matrix_file_contents.split(", ")

R_Matrix = []
for line in R_Matrix_file_contents:
    line = line[1:-1]
    line = list(line.split(" "))
    line = line[1:]
    
    line_float = [float(entry) for entry in line]
    R_Matrix.append(line_float)

R_Matrix = np.array(R_Matrix)


###-----Load Name References-----###
Name_References_file_name = sys.argv[8]
#R_Matrix_file_name = "/Ratio_Matrix/Ratio_Matrix_EFM_Nr_"+EFM_Nr+"_For.txt"
path = os.getcwd()+Name_References_file_name
Name_References_file = open(path,"r+")
Name_References_file_contents = Name_References_file.readlines()


max_tot_c = 0.5
nr_c_met = 0
for line in Name_References_file_contents:
    #print(line)
    if line[0] =="M":
        if "[e]" not in line:
            nr_c_met +=1
            if "h2o" in line:
                max_tot_c += 1

            if "biomass" in line:
                max_tot_c += 1

            if "PHB" in line:
                max_tot_c += 1







#########-----Algorithm------#########

# Constrain concentration ratios
# Use natural log
ratio_lim = np.log(np.array([
    [ 0.499, 50.1 ], # 0.5 < ATP / ADP < 50
    [ 0.00499, 0.501 ],  # 0.005 < NADH / NAD < 0.5
    [ 0.0499, 50.1 ], # 0.05 < NADPH / NADP < 50
    [ 0.099, 10.1 ]  # 0.1 < QH2 / Q < 10
]))


# Define function for random sampling of concentrations
def random_c(MetRange):
    sample = np.array([np.random.random() for n in range(0, MetRange.shape[0])])
    return sample * (MetRange[:,1] - MetRange[:,0]) + MetRange[:,0]

# Define function for checking if set is thermodynamically feasible
def df_ok(c,MDF):
    # Calculate delta G prime
    df = -(dG0 + RT * np.sum(np.transpose(S_Matrix) * c, 1))
    # Check if all driving forces higher than 0
    #print("df is:\n")
    #print(sum(df > MDF*0))
    # if not sum(df >= MDF*0.9) == df.shape[0]:
    #     print("It's the dGs!")
    return sum(df >= MDF*0.9) == df.shape[0]

# Define function for checking if set has acceptable ratios
def ratios_ok(c):
    #ratios = np.sum(R_Matrix.T * c, 1).reshape([ratio_lim.shape[0], 1])
    #print(ratios)
    ratios = np.sum(R_Matrix.T * c, 1).reshape([ratio_lim.shape[0], 1])
    min = np.sum(np.subtract(ratios, ratio_lim) >= 0, 0)[0] == ratios.shape[0]
    max = np.sum(np.subtract(ratios, ratio_lim) <= 0, 0)[1] == ratios.shape[0]
    # if not min or max:
    #     print("It's the ratios")
    return min and max

# Define function for checking that sum of concentrations is not too high (0.5 M)
def sum_ok(c, max_tot_c):
    #print("sum of all conc is:\n")
    #print(np.sum(np.exp(c)))
    



    ## Sum only intracellular metabolites
    return np.sum(np.exp(c_0[-nr_c_met:])) <= max_tot_c

# Define function that checks concentrations are within limits
def limits_ok(c):
    c_l = c.reshape([c.shape[0],1])
    min = np.sum(np.subtract(c_l, MetRange) >= 0, 0)[0] == c.shape[0]
    max = np.sum(np.subtract(c_l, MetRange) <= 0, 0)[1] == c.shape[0]
    # if not min or max:
    #     print("It's the ranges!")
    return min and max

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible(c,MDF,max_tot_c):
    return df_ok(c,MDF) and limits_ok(c) and ratios_ok(c) and sum_ok(c[2:],max_tot_c)
    print("Found feasible set!")

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible_final(c,MDF,max_tot_c):
    if not df_ok(c,MDF):
        print("It is the dG!")
    
    # if not ratios_ok(c):
    #     print("It is the ratios!")
    #     ratios = np.sum(R_Matrix.T * c, 1).reshape([ratio_lim.shape[0], 1])
    #     print(np.exp(ratios))

    if not limits_ok(c):
        print("It is the ranges!")

    return df_ok(c,MDF) and limits_ok(c) and ratios_ok(c) and sum_ok(c[2:],max_tot_c)
    print("Found feasible set!")

# Modify direction in order to get unstuck from concentration limits, a.k.a. The Unsticking Function TM
def unstick_direction(c, direction, MetRange):
    # Determine what metabolites are stuck at limits
    stuck = c.reshape((c.size,1)) == MetRange
    # Determine current signs of direction vector
    dirsign = np.sign(direction)
    # Pick a random sign for metabolites stuck at max
    max_sign = np.random.choice([-1,1], 1)
    # All directions for metabolites stuck at max must be the same sign
    dirsign[stuck[:,1] * dirsign != 0] = max_sign
    # All directions for metabolites stuck at min must be the opposite sign
    dirsign[stuck[:,0] * dirsign != 0] = -max_sign
    # Determine the directions that must change sign
    change_sign = dirsign != np.sign(direction)
    # Change the sign of directions that must change sign
    direction[change_sign] = direction[change_sign] * -1
    # Return the compatibility-modified "unstuck" direction vector
    return direction

# Define function for selecting a random direction
def random_direction(c):
    # Create a random vector of the same length as c
    direction = np.array([np.random.random() for n in range(0, c.shape[0])])
    # Subtract 0.5 to introduce negative directions
    direction = direction - 0.5
    # Set fixed concentration direction to zero
    direction[MetRange[:,1] - MetRange[:,0] == 0] = 0
    # Normalize length of direction vector
    normalized_direction = direction / np.linalg.norm(direction)
    return normalized_direction

# Define function to generate one feasible metabolite concentration set
def generate_feasible_c(MetRange, MDF,max_tot_c):
    c = random_c(MetRange) # Initialize c
    while not is_feasible(c, MDF,max_tot_c):
        c = random_c(MetRange) # Generate new c until feasible
    return c

# Determine minimum and maximum possible theta given concentration limits
def calculate_theta_hard_limit(c, direction, MetRange):
    # Find smallest fraction of direction that hits limit if added
    theta_max = np.vstack([
        (MetRange[:,1] - c)[direction != 0] / direction[direction != 0],
        (MetRange[:,0] - c)[direction != 0] / direction[direction != 0]
    ])
    #print(theta_max)
    theta_max = np.max(theta_max, 0)
    #print(theta_max)
    theta_max = min(theta_max[theta_max >= 0])
    #print(theta_max)
    # Find smallest fraction of direction that hits limit if subtracted
    theta_min = np.vstack([
        (c - MetRange[:,1])[direction != 0] / direction[direction != 0],
        (c - MetRange[:,0])[direction != 0] / direction[direction != 0]
    ])
    #print(theta_min)
    theta_min = np.max(theta_min, 0)
    #print(theta_min)
    theta_min = -min(theta_min[theta_min >= 0])
    #print(theta_min)
    return (theta_min, theta_max)

# Define function for determining minimum and maximum step length (theta)
def theta_range(c, direction, max_tot_c, precision=1e-3):
    # Define function for honing in on a theta limit
    def hone_theta(theta_outer, max_tot_c, theta_inner=0):
        if is_feasible(c + theta_outer * direction, MDF, max_tot_c):
            # If the outer theta is feasible, accept that solution
            theta_inner = theta_outer
        else:
            while abs(theta_outer - theta_inner) > precision:
                # Calculate a theta between outer and inner limits
                theta_cur = (theta_outer + theta_inner) / 2
                if is_feasible(c + theta_cur * direction, MDF, max_tot_c):
                    # Move outwards, set inner limit to current theta
                    theta_inner = theta_cur
                else:
                    # Move inwards, set outer limit to current theta
                    theta_outer = theta_cur
        # Return inner theta
        return theta_inner
    # Get hard limits on theta from concentrations
    theta_lim = calculate_theta_hard_limit(c, direction, MetRange)
    # Hone in on upper theta
    theta_upper = hone_theta(theta_lim[1],max_tot_c)
    # Hone in on lower theta
    theta_lower = hone_theta(theta_lim[0],max_tot_c)
    # Return results
    return [theta_lower, theta_upper]

# Define function for performing hit-and-run sampling within the solution space
def hit_and_run(S_Matrix, dG0, MetRange, ratio_lim, R_Matrix, n_samples, MDF, max_tot_c, precision=1e-3):
    # Generate starting point
    #c = generate_feasible_c(MetRange, MDF)
    #print("--- %s seconds to find the first feasible ---" % (time.time() - start_time))

    # Take starting point from Input
    c=c_0 

    # Set up concentration storage list
    fMCSs = [c]
    # Perform n steps
    for i in range(0, n_samples - 1):
        # Generate random direction
        direction = random_direction(c)
        # Make sure that the algorithm doesn't get stuck at the boundaries of the solution space
        direction_unstuck = unstick_direction(c, direction,MetRange)
        # Determine minimum and maximum step length
        theta = theta_range(c, direction_unstuck, max_tot_c, precision=precision)
        # Perform a random sampling of the step length
        theta = theta[0] + np.random.random() * (theta[1] - theta[0])
        # Perform step
        c = c + theta * direction
        # Ensure feasibility
        if not is_feasible_final(c,MDF,max_tot_c):
            print("Warning: Infeasible point reached.")
            break
        # Store concentration
        fMCSs.append(c)
    # Return list of concentrations
    return fMCSs

count = 0
final_out_Conc = ''
final_out_dG = ''
for c in hit_and_run(S_Matrix, dG0, MetRange, ratio_lim, R_Matrix, 5000, MDF, max_tot_c):
    # Print CSV in mM
    count+=1
    final_out_Conc = final_out_Conc + "fMCS"+str(count)+"," + ",".join([str(np.round(np.exp(x)*1000,3)) for x in c]) + "\n"
    df = -(dG0 + RT * np.sum(np.transpose(S_Matrix) * c, 1))
    final_out_dG = final_out_dG + "fMCS"+str(count)+"," + ",".join([str(np.round(df_1,3)) for df_1 in df]) + "\n"


Sampling_file_contents = final_out_dG

Sampling_file_contents = Sampling_file_contents.split('\n')
for line in Sampling_file_contents[:-1]:
    line_split = line.split(',')
    if line_split[0] == 'fMCS1':
        #print(line_split[0])
        line_split = line_split[1:]
        #print(line_split)
        line_split = [float(entry) for entry in line_split]
        Data_All = np.array(line_split)
    #elif line_split[0] != 'fMCS1':
    else:
        #print(line_split[0])
        line_split = line_split[1:]
        line_split = [float(entry) for entry in line_split]
        Data_fmc_Others = np.array(line_split)
        Data_All = np.vstack((Data_All,Data_fmc_Others))



## Calculate Median and MAD
medians = np.round(np.median(Data_All, axis=0),3)
final_out_Median = medians
#print(medians)
MADs = np.round(stats.median_abs_deviation(Data_All),3)
final_out_MAD = MADs
    #print(MADs)

    # Median_File = open("Medians_"+EFM_Nr+".txt","w")
    # np.savetxt(Median_File,medians)
    # Median_File.close()

    # MAD_File = open("MADs_"+EFM_Nr+".txt","w")
    # np.savetxt(MAD_File,MADs)
    # MAD_File.close()



# Output_File_Name_Conc = sys.argv[8]
# #Output_File = open("Sampling_Results_EFM_Nr_"+EFM_Nr+"_For_WT.txt","w")
# Output_File_Conc = open(Output_File_Name_Conc,"w")
# Output_File_Conc.write(final_out_Conc)
# Output_File_Conc.close()

# Output_File_Name_dG = sys.argv[9]
# Output_File_dG = open(Output_File_Name_dG,"w")
# Output_File_dG.write(final_out_dG)
# Output_File_dG.close()




Output_File_Name_Median = sys.argv[9]
#Output_File_Median = open(Output_File_Name_Median,"w")
path_Out_1 = os.getcwd()+Output_File_Name_Median
Output_File_Median = open(path_Out_1,"w")
np.savetxt(Output_File_Median,final_out_Median)
#Output_File_Median.write(final_out_Median)
Output_File_Median.close()

Output_File_Name_MAD = sys.argv[10]
#Output_File_MAD = open(Output_File_Name_MAD,"w")
path_Out_2 = os.getcwd()+Output_File_Name_MAD
Output_File_MAD = open(path_Out_2,"w")
np.savetxt(Output_File_MAD,final_out_MAD)
#Output_File_MAD.write(final_out_MAD)
Output_File_MAD.close()

