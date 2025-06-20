import numpy as np
from collections import Counter

import westpa
from westpa.analysis import Run

import matplotlib.pyplot as plt
from decimal import Decimal

#Jonathan Borowsky
#4/22/25

#---------------------------------------------------------------------------------------------------------
#load paired progress coordinates of every walker and its parent
#---------------------------------------------------------------------------------------------------------

#see https://westpa.readthedocs.io/en/latest/documentation/analysis/index.html
#and https://westpa.readthedocs.io/en/latest/users_guide/hdf5.html

def h5_2_pcs(h5path, miniter, maxiter):
    
    run = Run.open(h5path)

    #set maximum iteration automatically
    if maxiter == -1:
        maxiter = run.num_iterations

    #collect matching arrays of the progress coordinate of each walker and its parent in the previous round
    #i.e. parent_pcs[i] contains the pc of the parent of the walker whose pc is recorded at pcs[i]
    pcs = []
    parent_pcs = []
    weights = []     #for parent walkers

    #initial progress coordinate (only collected if miniter <= 1)
    pc_init = []

    #loop over all WESTPA rounds
    for iteration in run:

        #check if current iteration is in the correct range
        if iteration.number < miniter:
            if iteration.number == 1:
                #initial pc
                pc_init = [[wa.parent.pcoord] for wa in iteration.walkers][0][-1]            
            continue
        elif iteration.number >= maxiter:
            break

        #collect PCs for current round
        pcs.append(iteration.pcoords)

        #collect PCs for parents of walkers in current round
        #note that these are collected as lists of 2d arrays rather than 3d arrays
        #    they are stacked rather than concatenated at the end to achieve the same output shape
        if iteration.number == 1:
            #for the first round there is no parent walker in the regular sense 
            #     but the pcoord is recorded under a slightly different variable name
            #     weight is similarly recorded differently
            parent_pcs.append(np.stack([[wa.parent.pcoord, wa.parent.pcoord] for wa in iteration.walkers]))
            weights.append([wa.parent.basis_state.probability for wa in iteration.walkers])
            
            #initial pc
            pc_init = parent_pcs[-1][-1]
            
        else:
            #not everything appears in .__dict__
            parent_pcs.append(np.stack([wa.parent.pcoords for wa in iteration.walkers]))
            weights.append([wa.parent.weight for wa in iteration.walkers])


    return pcs, parent_pcs, weights, pc_init, maxiter


#---------------------------------------------------------------------------------------------------------
#digitize and flatten a list of walker PCs
#---------------------------------------------------------------------------------------------------------

#returns a list of the bin assignments of the coordinates in pcs
def digitize_flatten_trj(pcs, binbounds, n_discrete_pc_vals):

    #put walkers in discrete bins
    connectivity_bins = np.digitize(pcs[:,-1,0], binbounds)
    #path_bins = [int(round(i)) for i in pcs[:,-1,1]]

    #maximum number of bins along the continuous PC axis; these may not all be occupied
    #n_connectivity_bins = len(binbounds)+1

    #flatten data in 2d PC space into a 1d PC space consisting of blocks along the discrete PC
    #bins_flattened = [pb*n_connectivity_bins + cb for pb, cb in zip(path_bins, connectivity_bins)]

    return connectivity_bins #bins_flattened


#---------------------------------------------------------------------------------------------------------
#digitize and flatten parent and child walker PCs and return a list of discrete transitions
#---------------------------------------------------------------------------------------------------------

def pcs_2_transitions(pcs, parent_pcs, binbounds, n_discrete_pc_vals):

    walker_bins_flattened = digitize_flatten_trj(pcs, binbounds, n_discrete_pc_vals)
    parent_bins_flattened = digitize_flatten_trj(parent_pcs, binbounds, n_discrete_pc_vals)
    
    return np.stack((parent_bins_flattened, walker_bins_flattened)).transpose()


#---------------------------------------------------------------------------------------------------------
#get a list of transitions and their weights using the methods above
#---------------------------------------------------------------------------------------------------------

def h5_2_transitions(h5path, miniter, maxiter, binbounds, n_discrete_pc_vals):

    pcs, parent_pcs, weights, pc_init, maxiter = h5_2_pcs(h5path, miniter, maxiter)
    transitions = pcs_2_transitions(np.concatenate(pcs), np.concatenate(parent_pcs), binbounds, n_discrete_pc_vals)

    return transitions, weights, pc_init, maxiter


#---------------------------------------------------------------------------------------------------------
#weight walkers based on the MSM weights of the bins they occupy and the number of walkers in each bin 
#so that the total weight of all walkers in each bin equals that bin's MSM probability
#---------------------------------------------------------------------------------------------------------

def get_walker_mi_weights(pcs, binbounds, n_discrete_pc_vals, pyem):

    #put walkers in discrete bins
    connectivity_bins = np.digitize(pcs[:,0,0], binbounds)
    path_bins = [int(round(i)) for i in pcs[:,0,1]]

    #maximum number of bins along the continuous PC axis; these may not all be occupied
    n_connectivity_bins = len(binbounds)+1

    #flatten data in 2d PC space into a 1d PC space consisting of blocks along the discrete PC
    bins_flattened = [pb*n_connectivity_bins + cb for pb, cb in zip(path_bins, connectivity_bins)]

    bin_counts = Counter(bins_flattened)

    # print(np.where(pyem.active_set == 231)[0])
    # print(np.where(pyem.active_set == 231)[0][0])
    try:
        walker_mi_weights = [pyem.stationary_distribution[np.where(pyem.active_set == walker_bin)[0][0]]/bin_counts[walker_bin] for walker_bin in bins_flattened]
        
    except IndexError as e:
        print("The following error probably arises from insufficient sampling. This leads to some states for which water wire and/or solvation data are present being ergodically trimmed from the MSM, causing an index mismatch. All walkers will be given equal weights instead of MSM-derived ones:\n")
        print(e)
        print("")
        walker_mi_weights = [1/(len(bins_flattened)*bin_counts[walker_bin]) for walker_bin in bins_flattened]
        
    return walker_mi_weights


#--------------------------------------------------------------
# Plot PyEMMA markov state model energies
#--------------------------------------------------------------
#TODO write method spec.

# Map PyEMMA MSM state indices back to progress coordinate values, convert them to energies,
# and plot a line for the energy as a function of the continuous PC at each discrete PC value

# set std=True for bayesian MSM only
def plot_2d_pc_webins(pyem, nbins, binset, pclims, discrete_pc_vals, pcinit, threshold, std=False, plottitle="", zeromin=False):

    bins = binset

    # -------------------------------------------------------
    # convert populations to energies
    energy = [-np.log(p) for p in pyem.stationary_distribution]
    min_energy = min(energy)
    energy = [ei - min_energy for ei in energy]

    # -------------------------------------------------------
    # align water wire PC values for wires through the three subunits and plot them against each other
    # 'align' here means to account for the index shift
    # due to PyEMMA discarding disconnected energy landscape components,
    # causing it to output fewer equilibrium probabilities than the original number of bins

    # empty list of bin indices for each PC 2 value
    # note that these are the bins used by h5_2_transitions_forpyemma() to discretize the PC for MSM construction
    # which need not match the original westpa bins
    bin_vals = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state energies and standard deviations thereof for each PC 2 value
    energies = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state occupancies and standard deviations thereof for each PC 2 value
    occupancies = [[] for x in range(discrete_pc_vals)]

    # pyem.active_set is a list of the input indices of the bins which are part of the largest connected component
    # state i in the MSM returned by PyEMMA therefore corresponds to bin pyem.active_set[i]
    # note that these are bin indices in the flattened data sent to PyEMMA, not just PC 1 bin indices
    # also note that they are 1-indexed
    for i, b in enumerate(pyem.active_set):
        # the value of PC 2
        #  Though naturally distributed in 2D PC1-PC2 space, bins were flattened
        #  for transition matrix (and MSM) construction, so the transition matrix provided to PyEMMA
        #  has dimensions of (discrete_pc_vals*nbins) x (discrete_pc_vals*nbins),
        #  with all the bins sharing each PC2 value grouped together in order.
        #  The first nbins PyEMMA input states therefore have the lowest PC 2 value,
        #  the next nbins PyEMMA input states the second lowest, etc.
        #  I think 1 is subtracted because pyem.active_set treats the input states as 1-indexed.
        #  This assumes a 0-indexed integer PC 2 value.
        pc2_val = (b - 1) // nbins

        # separate MSM data by PC 2 value

        # determine the bin index of the current bin along PC 1 and add it to the list
        # note that (b-1)%nbins = b-1-nbins*pc2_val
        bin_vals[pc2_val].append((b - 1) % nbins)

        # add MSM state energies, occupancies, and standard deviations thereof to lists for subsequent plotting
        energies[pc2_val].append(energy[i])
        #if std: std_energies[pc2_val].append(std_energy[i])

        occupancies[pc2_val].append(pyem.stationary_distribution[i])
        #if std: std_occupancies[pc2_val].append(pyem.sample_std('stationary_distribution')[i])

    # -------------------------------------------------------
    # Plot energies as a function of the progress coordinate
    # specifically plot energies as a function of PC 1 for each value of PC 2
    # also calculate total occupancy below the PC 1 threshold for each value of PC 2

    bin_x_all = []
    
    q_below_total = 0
    
    for i in range(discrete_pc_vals):
                
        bin_x = [bins[b] for b in bin_vals[i]]
        #the bin at 1.0 has the combined weights of everything past PC0=1
        # plt.plot(bin_x, energies[i]) 
        bin_x_all.append(bin_x)
        # # occupancy below threshold, which should be insensitive to bin size
        # print(f"path {i}")

        # # total occupancy of bins below threshold
        # q_below = sum([o for o, x in zip(occupancies[i], bin_x) if x < threshold])
        # print(f"{'%.3E' % Decimal(q_below)} occupancy below {threshold} nm")
        # q_below_total += q_below

        # standard deviation of total occupancy of bins below threshold
        # variances add, standard deviations do not
        # TODO we should verify whether our assumption of independent bin variance is correct
        # if std:
        #     var_below = sum([o ** 2 for o, x in zip(std_occupancies[i], bin_x) if x < threshold])
        #     print(f"{'%.1E' % Decimal(np.sqrt(var_below))} standard deviation in occupancy below {threshold} nm")

    # -------------------------------------------------------
    # plot the threshold PC 1 value and the starting structure PC
    # and set plot labels and legends and such

    # # matplotlib already selects reasonable y axis bounds
    # ylims = plt.gca().get_ylim()

    # # get the standard matplotlib colors in order
    # prop_cycle = plt.rcParams['axes.prop_cycle']
    # colors = prop_cycle.by_key()['color']

    # # plot the progress coordinate of the starting structure
    # plt.plot([pcinit[0], pcinit[0]], [ylims[0], ylims[1]], color=colors[round(pcinit[0])], linestyle='dotted')

    # plt.legend([f"path {i + 1}" for i in range(discrete_pc_vals)] + ["starting structure"])

    # # plot a line denoting the PC 1 threshold
    # plt.plot([threshold, threshold], [ylims[0], ylims[1]], color='black', linestyle='dashed')


    # #binwidth = .5
    # #plt.xlim(pclims[0]-binwidth/2, pclims[1]-binwidth/2)
    # #plt.xlim(0.25, 1)
    # # plotting lines alters the ylims, which is not desired
    # #plt.ylim(ylims[0], ylims[1])
    # #plt.ylim(0, 12)
    # plt.xlabel("maximum inter-oxygen distance along water path (nm)")
    # plt.ylabel("free energy (kT)")

    # if plottitle != "":
    #     plt.savefig(plottitle, dpi=600)
    #     print("saved")
    #     #plt.savefig(plottitle+".svg", dpi=600, format="svg")

    # plt.show()
    
    # print("""\nWARNING
    # Beware that as water wires of similar lengths need not be adjacent 
    # or even remotely near each other in configuration space, 
    # the height of the energy landscape between conducting and experimental states need not 
    # accurately represent the energy barrier to the formation of those conducting states. 
    # haMSMs will be used to more accurately determine the relative kinetics of formation 
    # of different conducting states in future work.""")

    return q_below_total, bin_x_all, energies
