#Jonathan Borowsky
#02/06/25
#Grabe lab

import numpy as np
import h5py

#------------------------------------------------------------------------------------------------------------
#identify the ancestors of a westpa walker
#------------------------------------------------------------------------------------------------------------

#copied from westpa-barrier-crossing-v3 by way of other notebooks
#TODO: this is used in many places and needs to be turned into a package or importable .py file

def walker_ancestors(h5path, walker_round, walker_num):
    
    #load h5 file
    with h5py.File(h5path, 'r') as f:
        
        pcoords = []
        walker_ids = []
        
        #determine names of westpa rounds
        #note that round 1's name is at index 0
        iterations = [iter for iter in f["iterations"]]

        current_round = walker_round
        
        for iter_ind in range(walker_round):
        
            #this extracts only the iteration name, not the iteration data
            iter_name = iterations[current_round-1] 
            #using the iteration name to extract the data
            
            #log walker ID and progress coordinate
            #zeros are for trimming excess nested array layers
            pcoords.append(f["iterations"][iter_name]["pcoord"][walker_num][0][0])
            walker_ids.append(walker_num)
            
            #update round and walker numbers
            current_round -= 1
            if current_round == 0:
                break
            walker_num = f["iterations"][iter_name]["seg_index"][walker_num][1]

    walker_ids.reverse()
    pcoords.reverse()
    return [walker_ids, pcoords]


#------------------------------------------------------------------------------------------------------------
#             convert microscopic free energy profiles to macroscopic binding coefficients 
#             this version is a test for unbound states in isotropic 3d solution
#------------------------------------------------------------------------------------------------------------

#x = coordinates of states
#fe = free energies of states from simulation, in k_B*T
#macrostate_classifier: function that takes a state coordinate from x and returns 0 if the state is bound and a different number otherwise
#T: temperature in Kelvin
#V_solvent: volume of solvent in the simulation in liters
#    For a single solvent phase this should be V_box - V_protein, which will be approximately V_box in most cases. 
#    For an NPT ensemble the average box volume should be used.
#    It is assumed that the volume is approximately the same in bound and unbound states.
#reference_concentration: a reference concentration in Molar (mol/liter) 
#    Used for converting from macroscopic binding equilibrium coefficients to binding free energies.

def fe_to_kd_3d_solution(x, fe, macrostate_classifier, T, V_solvent, reference_concentration=1):
    #constants
    k_B = 1.380649*10**-23   #J/K      #Boltzmann's constant
    N_A = 6.02214076*10**23  #mol**-1  #Avogradro's number

    #partition function
    Z = sum([np.exp(-fei) for fei in fe])

    #bound fraction
    f_b = (1/Z)*sum([np.exp(-fei) for x, fei in zip(x, fe) if macrostate_classifier(x) == 0])

    #macroscopic binding and dissociation coefficients
    K_bind = N_A*V_solvent*f_b/(1-f_b)**2  #L/mol
    K_d = 1/K_bind                         #mol/L = Molar
    
    #macroscopic binding free energy per mole
    G_b = -N_A*k_B*T*np.log(reference_concentration*K_bind)

    return K_d, G_b


#------------------------------------------------------------------------------------------------------------
#classify bound and unbound states; for use in free energy to binding coefficient conversion
#------------------------------------------------------------------------------------------------------------

def macrostate_classifier(x):
    if x < 5:
        return 0
    else:
        return 1

