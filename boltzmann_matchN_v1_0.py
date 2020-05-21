###################################################################################
###################################################################################
###                                                                             ###
###     Simulation of Spontaneous Emergence of Boltzmann Distribution           ###
###                             Version 1.0                                     ###
###                                                                             ###
###              by Gregory Jameson and Rafael Bruschweiler                     ###
###                                                                             ###
###                     The Ohio State University                               ###
###                                                                             ###
### Simulates the spontaneous emergence of a Boltzmann-like distribution in     ###
### a system consisting of N distinguishable molecules through: 	               ###
###       - Equal a priori probabilities of microstates                         ###
###       - Conservation of total energy                                        ###
###       - This version of the program tries to match energy of                ### 
###         Nth molecule so that total energy is conserved                      ###
###         Note: All energies are in units of epsilon                          ###
###                                                                             ###
### For more explanations see:                                                  ###
###     G. Jameson and R. Bruschweiler (submitted for publication)              ###
###                                                                             ###
### This software is provided "as is" and any expressed or implied              ###
### warranties are disclaimed.                                                  ###
###################################################################################
###################################################################################

import numpy as np
import matplotlib.pyplot as plt

Trials = 1000000 # Set number of trials to run; start with small numbers to reduce compute time and then gradually increase

N = 15;  # Number of molecules 
E_total = 40;  # Total energy of all molecules combined in units of epsilon
E_min = 0; # Minimal energy of a molecule
E_max = 10; # Maximal energy of a molecule

# Optional: Insert HERE line for the examples at the bottom: N = ...; 
E_levels = np.arange(E_min,E_max+1) # Energy levels that can be occupied by an individual molecule (including degeneracies)

# Initialize loop to generate energy configuration where each molecule is assigned a random energy from vector E_levels 
count = 0;
stat = np.zeros(E_levels.size);

# Start energy configuration generation loop and select acceptable configurations
for k in range(Trials):
    E_select = np.random.randint(E_levels.size,size=(N,)) # Generate energy configuration by assigning random energy to each molecule 
    E_missing = E_total - sum(E_levels[E_select]) # Determine the error in total energy
    E_select_fill = np.argwhere(E_missing==E_levels) # Generate the energy in the final molecule if it matches an energy of vecotr E_levels
    
    # If the missing energy is in the set of available energy levels then accept configuration, otherwise not
    if E_select_fill.size == 1:
        E_select = np.concatenate( (E_select, E_select_fill[0]) ) # Add missing energy level to set of energy levels
        for E_this in E_select:
            stat[E_this] = stat[E_this] + 1 # Update the occurrence of the energy levels of accepted configurations
        count = count + 1  # Update the count of number of successful trials

pop = stat/sum(stat)  # Convert number of energy occurrences to normalized energy level populations

# The following evaluations of energy and energy fluctuations may be useful for further analysis
E_tot_avg = sum(np.multiply(stat,E_levels))/count  # Average total energy of all accepted configurations
E_mol_avg = sum(np.multiply(pop,E_levels))  # Average energy per molecule
E2_levels = np.power(E_levels,2)  # Needed to compute standard deviation
E2_mol_avg = sum(np.multiply(pop,E2_levels))  # Average energy per molecule
E_mol_std = np.sqrt(E2_mol_avg - E_mol_avg**2)  # Standard deviation of energy per molecule


# Display final results
print('Number of successful trials: ')
print(count) # Display number of successful trials

# Plot energy distribution of individual molecules
plt.clf()
plt.rcParams.update({'font.size':16})
plt.plot(E_levels,pop,'o',markersize=10,markerfacecolor='b',linewidth=3)
plt.title('Molecular energy distribution')
plt.xlabel('Energy of a molecule in units of epsilon')
plt.ylabel('Probability p_i')
plt.show()

# Plot energy distribution as a semilog plot to assess whether thermodynamic limit has been approximately reached 
plt.clf()
plt.rcParams.update({'font.size':16})
plt.semilogy(E_levels,pop,'o',markersize=10,markerfacecolor='b',linewidth=3)
plt.title('Molecular energy distribution')
plt.xlabel('Energy of a molecule in units of epsilon')
plt.ylabel('Probability p_i')
plt.show()

'''
The following input parameters are suggestions to try out.   
Insert lines below "N = ...;" before E_levels = [E_min:E_max]; 
Important: try other inputs on your own and analyze the results

For Figure 1 of the paper:
N = 4; E_total = 6; E_error = 0; E_min = 0; E_max = 6; Trials = 1000000;

For Figure 2 of the paper:
N = 30; E_total = 900; E_error = 0; E_min = 0; E_max = 100; Trials = 100000000;

Other suggestion:
N = 10; E_total = 12; E_error = 0; E_min = 0; E_max = 6; Trials = 10000000;

Other suggestion:
N = 15; E_total = 40; E_error = 1; E_min = 0; E_max = 10; Trials = 10000000;

Other suggestion:
N = 4; E_total = 18; E_error = 0; E_min = 0; E_max = 6; Trials = 1000000;
Try to rationalize the result. What real world systems can show such behavior? 
'''