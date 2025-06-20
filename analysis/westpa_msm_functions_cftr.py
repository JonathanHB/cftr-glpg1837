#Jonathan Borowsky
#02/06/25
#Grabe lab

import os
import sys

import numpy as np
import h5py
import pyemma
import time

import matplotlib.pyplot as plt
from decimal import Decimal

#--------------------------------------------------------------
# Extract and discretize trajectory segments for pyemma from the westpa west.h5 file
#--------------------------------------------------------------

#TODO write method specification here

def recursive_npyloader(files, npypath, x):
    
    for j in range(x, 0, -10):
        if f"pc_data_{j}_v1.npy" in files:
            try:
                print(f"loading {npypath}/pc_data_{j}_v1.npy")
                data_seg = np.load(f"{npypath}/pc_data_{j}_v1.npy")
                
                #print(min(data_seg[:,0]))
                if min(data_seg[:,0])>1:
                    data_loaded = recursive_npyloader(files, npypath, int(min(data_seg[:,0])-1))
                else:
                    data_loaded = []
                
                data_loaded.append(data_seg)
                
                break
            except Exception as e:
                print(f"skipped {npypath}/pc_data_{j}_v1.npy with exception {e}")

    return data_loaded

    
def npyloader(npypath):

    files = os.listdir(npypath)
    #print(files)
    rf_data = recursive_npyloader(files, npypath, 10000)
    rf_data = np.concatenate(rf_data)
    
    return rf_data


def h5_2_transitions_forpyemma_webins(h5path, npypath, minround=0, maxround=-1, allcomplete=False, discrete_pc_vals=3,
                               binset = [-1,0,1], n_walkers=6):

    minsofar = 999
    minpcwalker = []
    minpcround = []
    
    pcs_all0 = []
    parent_pcs_all0 = []  # for calculating histogram/state limits

    pcs_all1 = []
    parent_pcs_all1 = []  # for calculating histogram/state limits

    init_pc_01 = [-1, -1]  # PC of the single seed structure; these are dummy values for cases where minround != 0

    use_rf_pc = False
    if npypath != "":
        rf_data = npyloader(npypath)
        use_rf_pc = True

        #print(rf_data)
        # print(np.max(rf_data, axis=0))
    
    # load h5 file
    with h5py.File(h5path, 'r') as f:

        # determine number of westpa rounds and set number to use accordingly
        iterations = [iter for iter in f["iterations"]]
        maxiter = len(iterations)
        if maxround == -1:
            if allcomplete:
                maxround = maxiter
            else:  # if the last round is not complete, omit it as data have yet to be written
                maxround = maxiter - 1

        if use_rf_pc:
            maxround = int(min([maxround, rf_data[-1][0]]))
            # print(maxround)
        
        print(f"loading data for {maxround - minround} westpa rounds")

        # used to calculate transitions, saved to avoid having to extract these twice for every iteration
        lastiter_pcs = []

        for iter_ind in range(minround, maxround):
            #print(iter_ind)
            # this extracts only the iteration name, not the iteration data
            iter_name = iterations[iter_ind]
            # using the iteration name to extract the data
            iter_data = f["iterations"][iter_name]

            # extract progress coordinate values and weights from array structure
            if not use_rf_pc:
                pcs0 = [i[-1][0] for i in iter_data["pcoord"]]
            else:
                pcs0 = [rfd[2] for rfd in rf_data if rfd[0] == iter_ind+1]
                # if len(pcs0) > 0:
                #     if min(pcs0) < minsofar and min(pcs0) >= 0:
                #         minsofar = min(pcs0)
                #         minpcwalker = np.argmin(pcs0)
                #         minpcround = iter_ind+1
                #print(len(pcs0))
                #print([rfd for rfd in rf_data if rfd[0] == iter_ind+1])
                
            pcs1 = [0 for i in iter_data["pcoord"]] #1DMOD
            #[i[0][1] for i in iter_data["pcoord"]]

            # when not starting from the beginning of the westpa walker tree,
            # get the progress coordinates of the first round used
            # and then proceed to the next round (as a single round has no transitions)
            if iter_ind == minround and minround != 0:
                lastiter_pcs0 = pcs0
                lastiter_pcs1 = pcs1
                continue

            skipround = False
            # get walker parent progress coordinates
            if iter_ind != 0:  # for the first round parent ids are negative
                # get walker parent IDs
                parent_inds = [i[1] for i in iter_data["seg_index"]]
                #to omit a corrupted round
                if max(parent_inds) > len(lastiter_pcs0): #iter_ind == 494 and 
                    #print(f"skipping iteration {iter_ind}, missing data")
                    skipround = True
                else:
                    parent_pcs0 = [lastiter_pcs0[j] for j in parent_inds]
                    parent_pcs1 = [lastiter_pcs1[j] for j in parent_inds]

            else:
                # starting structure progress coordinate
                if not use_rf_pc:
                    init_pc0 = [i for i in iter_data["ibstates"]["bstate_pcoord"]][0][0]
                else:
                    init_pc0 = rf_data[0,2]
                    
                parent_pcs0 = [init_pc0 for i in range(n_walkers)]

                init_pc1 = 0 #[i for i in iter_data["ibstates"]["bstate_pcoord"]][0][1] #1DMOD
                parent_pcs1 = [init_pc1 for i in range(n_walkers)]

                init_pc_01 = [init_pc0, init_pc1]

            #if len(parent_pcs0) != len(pcs0):
                #print(f"skipping iteration {iter_ind}; mismatched data")
                # #print(parent_pcs0)
                #print(len(parent_pcs0))
                # #print(pcs0)
                #print(len(pcs0))
                # print("---------------------")
                
                # if iter_ind == 3:
                # import sys
                # sys.exit(0)

            # assemble arrays of progress coordinates shifted by 1 round
            if (not skipround) and (len(parent_pcs0) == len(pcs0)) and (-1.0 not in parent_pcs0) and (-1.0 not in pcs0):
                pcs_all0 += pcs0
                parent_pcs_all0 += parent_pcs0
    
                pcs_all1 += pcs1
                parent_pcs_all1 += parent_pcs1

            # update last round's progress coordinates
            lastiter_pcs0 = pcs0
            lastiter_pcs1 = pcs1

    #print(pcs_all0)
    #print(parent_pcs_all0)
    
    # progress coordinate 0 overall range
    endpad = 10 ** -6
    pcmin = min(pcs_all0 + parent_pcs_all0) - endpad
    pcmax = max(pcs_all0 + parent_pcs_all0) + endpad

    trjs_binned_all = []

    trjs_binned_all0 = []

    #for nbins in binrangeobj:
    # define bin edges
    # This could be done with one arange command but I'm afraid that
    # rounding/floating point issues could cause the last value to go missing
    #intrange = np.arange(nbins * discrete_pc_vals + 1)
    #bins = intrange / nbins * (pcmax - pcmin) + pcmin
    # print(len(bins))
    
    #bin_offset_pc1 = 100 #must be much larger than the pc0 range and maximum pc0 values
    
    bin_range = max(binset) - min(binset)
    # print("------------------")
    # print(min(pcs_all0 + parent_pcs_all0))
    # print(max(pcs_all0 + parent_pcs_all0))
    #pc0range = pcmax - pcmin
    
    bins = []
    for i in range(discrete_pc_vals):
        bins += [bb + bin_range*i for bb in binset]

    # print(bins)
    # separate the progress coordinates by pc1 (which pseudotrimeric unit they go through)
    pcs_all_1d = [pc0 + np.round(pc1) * bin_range for pc0, pc1 in zip(pcs_all0, pcs_all1)]
    parent_pcs_all_1d = [pc0 + np.round(pc1) * bin_range for pc0, pc1 in zip(parent_pcs_all0, parent_pcs_all1)]

    #plt.hist2d(pcs_all_1d, parent_pcs_all_1d, bins=(40,40))
    #plt.show()
    
    # stack paired parent and child pcs
    pcs_all_parent_pcs_all = np.stack((parent_pcs_all_1d, pcs_all_1d))

    # bin pcs into states
    trjs_binned = np.digitize(pcs_all_parent_pcs_all, bins).transpose()

    trjs_binned_all.append(trjs_binned)

    #collect transitions along PC0 regardless of PC1 value
    pcs_all_parent_pcs_all0 = np.stack((parent_pcs_all0, pcs_all0))
    trjs_binned0 = np.digitize(pcs_all_parent_pcs_all0, bins).transpose()

    trjs_binned_all0.append(trjs_binned0)

    # print(minsofar)
    # print(minpcround)
    # print(minpcwalker)

    return trjs_binned_all, (min(binset), max(binset)), init_pc_01, trjs_binned_all0, (pcmin, pcmax)



#--------------------------------------------------------------
# Plot PyEMMA markov state model energies
#--------------------------------------------------------------
#TODO write method spec.

# Map PyEMMA MSM state indices back to progress coordinate values, convert them to energies,
# and plot a line for the energy as a function of the continuous PC at each discrete PC value

# set std=True for bayesian MSM only
def plot_2d_pc_webins(pyem, binset, pclims, discrete_pc_vals, pcinit, threshold, std=False, plottitle=""):
    # -------------------------------------------------------
    # generate bin boundary coordinates along PC 1 (aka PC 0 if 0 indexing)
    # a number of integers equal to the maximum total number of bins
    #intrange = np.arange(discrete_pc_vals * nbins)
    # width of bins along PC 1
    #binwidth = (pclims[1] - pclims[0]) / nbins
    # bin boundary coordinates along PC 1 extending for thrice the actual PC 1 range
    bins = binset #intrange * binwidth + pclims[0]

    # -------------------------------------------------------
    # convert populations to energies
    energy = [-np.log(p) for p in pyem.stationary_distribution]

    # TODO figure out why this throws an error (possibly the lower ends of error bars being below 0) and fix it
    # TODO I'm not certain this is the right way to convert standard deviations from populations to energies
    if std: std_energy = [-np.log(p - stdp) + np.log(p + stdp) for p, stdp in
                          zip(pyem.stationary_distribution, pyem.sample_std('stationary_distribution'))]

    # -------------------------------------------------------
    # align water wire PC values for wires through the three subunits and plot them against each other
    # 'align' here means to account for the index shift
    # due to PyEMMA discarding disconnected energy landscape components,
    # causing it to output fewer equilibrium probabilities than the original number of bins
    # no attempt to make the energy values themselves match is involved <--?

    # empty list of bin indices for each PC 2 value
    # note that these are the bins used by h5_2_transitions_forpyemma() to discretize the PC for MSM construction
    # which need not match the original westpa bins
    bin_vals = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state energies and standard deviations thereof for each PC 2 value
    energies = [[] for x in range(discrete_pc_vals)]
    if std: std_energies = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state occupancies and standard deviations thereof for each PC 2 value
    occupancies = [[] for x in range(discrete_pc_vals)]
    if std: std_occupancies = [[] for x in range(discrete_pc_vals)]

    # pyem.active_set is a list of the input indices of the bins which are part of the largest connected component
    # state i in the MSM returned by PyEMMA therefore corresponds to bin pyem.active_set[i]
    # note that these are bin indices in the flattened data sent to PyEMMA, not just PC 1 bin indices
    # also note that they are 1-indexed <-- nevermind a zero has been sighted; seemingly they are 0-indexed
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
        pc2_val = (b - 1) // len(binset)

        # separate MSM data by PC 2 value

        # determine the bin index of the current bin along PC 1 and add it to the list
        # note that (b-1)%nbins = b-1-nbins*pc2_val
        bin_vals[pc2_val].append((b - 1) % len(binset))

        # add MSM state energies, occupancies, and standard deviations thereof to lists for subsequent plotting
        energies[pc2_val].append(energy[i])
        if std: std_energies[pc2_val].append(std_energy[i])

        occupancies[pc2_val].append(pyem.stationary_distribution[i])
        if std: std_occupancies[pc2_val].append(pyem.sample_std('stationary_distribution')[i])

    #print(binset)
    
    # -------------------------------------------------------
    # Plot energies as a function of the progress coordinate
    # specifically plot energies as a function of PC 1 for each value of PC 2
    # also calculate total occupancy below the PC 1 threshold for each value of PC 2

    bin_x = [bins[b] for b in bin_vals[0]]
    #print(bin_x)
    return [bin_x, energies[0]]

    
    # return_data = []
    
    # for i in range(discrete_pc_vals):
                
    #     bin_x = [bins[b] for b in bin_vals[i]]

    #     # plt.errorbar(bin_x, occupancies[i], yerr=std_occupancies[i])
    #     if std:
    #         plt.errorbar(bin_x, energies[i], yerr=std_energies[i])
    #     else:
    #         #plt.plot(bin_x, energies[i]) #the bin at 1.0 has the combined weights of everything past PC0=1

    #         return_data = [bin_x, energies[i]]
        
    #     # occupancy below threshold, which should be insensitive to bin size
    #     print(f"path {i}")

    #     # total occupancy of bins below threshold
    #     q_below = sum([o for o, x in zip(occupancies[i], bin_x) if x < threshold])
    #     print(f"{'%.3E' % Decimal(q_below)} occupancy below {threshold} nm")

    #     # standard deviation of total occupancy of bins below threshold
    #     # variances add, standard deviations do not
    #     # TODO we should verify whether our assumption of independent bin variance is correct
    #     if std:
    #         var_below = sum([o ** 2 for o, x in zip(std_occupancies[i], bin_x) if x < threshold])
    #         print(f"{'%.1E' % Decimal(np.sqrt(var_below))} standard deviation in occupancy below {threshold} nm")

    # # -------------------------------------------------------
    # # plot the threshold PC 1 value and the starting structure PC
    # # and set plot labels and legends and such

    # # matplotlib already selects reasonable y axis bounds
    # ylims = plt.gca().get_ylim()

    # # get the standard matplotlib colors in order
    # # prop_cycle = plt.rcParams['axes.prop_cycle']
    # # colors = prop_cycle.by_key()['color']

    # # plot the progress coordinate of the starting structure
    # #plt.plot([pcinit[0], pcinit[0]], [ylims[0], ylims[1]], color=colors[round(pcinit[1])], linestyle='dotted')

    # #plt.legend([f"path {i + 1}" for i in range(discrete_pc_vals)] + ["starting structure"])

    # # plot a line denoting the PC 1 threshold
    # #plt.plot([threshold, threshold], [ylims[0], ylims[1]], color='black', linestyle='dashed')

    # binwidth=1
    # plt.xlim(pclims[0]-binwidth/2, pclims[1]+binwidth/2)
    # #plt.xlim(0.25, 1)
    # # plotting lines alters the ylims, which is not desired
    # plt.ylim(ylims[0], ylims[1])
    # #plt.xlabel("maximum inter-oxygen distance along water path (nm)")
    # #plt.ylabel("free energy (kT)")

    # if plottitle != "":
    #     plt.savefig(plottitle+".svg", dpi=600, format="svg")

    # print("""\nBeware that as water wires of similar lengths need not be adjacent 
    # or even remotely near each other in configuration space, 
    # the height of the energy landscape between conducting and experimental states need not 
    # accurately represent the energy barrier to the formation of those conducting states. 
    # haMSMs will be used to more accurately determine the relative kinetics of formation 
    # of different conducting states in future work.""")

    # return return_data



#--------------------------------------------------------------
# Construct a bayesian markov state model from a westpa h5 file and plot the energies
# runs the two methods above
#--------------------------------------------------------------
#TODO write method spec.

def build_pyemma_msm_webins(h5path, npypath, minround, maxround, n_discrete_pc_vals, binrangeobj, threshold, n_walkers=6, plot_pyemma_bayesian_error_bars=True, savefigname=""):

    t1 = time.time()
    # get all transitions, with a number of equally spaced PC1 bins as specified by binrangeobj
    trjs, pclims, pcinit, trjs_binned_all0, pcextremes = h5_2_transitions_forpyemma_webins(h5path, npypath, minround, maxround, discrete_pc_vals=n_discrete_pc_vals, binset=binrangeobj, n_walkers=n_walkers)
    t2 = time.time()
    print(f"loaded data for {trjs[0].shape[0]} transitions in {t2-t1} seconds")

    # build msm for each number of bins

    #nbins = [len(binrangeobj)]

    trj = trjs[0]
    #print(trj)
    #for trj, nb in zip(trjs, nbins):
        # build MSM
        # note that adding a semicolon does not suppress Intel MKL warnings here
    t3 = time.time()
    pyem = pyemma.msm.estimate_markov_model(list(trj), lag=1, reversible=True)  # , mincount_connectivity=0
    t4 = time.time()
    print(f"built msm in {t4-t3} seconds")
    
    # plot MSM energies
    plot_data = plot_2d_pc_webins(pyem, binrangeobj, pcextremes, n_discrete_pc_vals, pcinit, threshold, plot_pyemma_bayesian_error_bars, savefigname)

    # plot transition matrix (this looks like garbage for the 2d PCs since it's been flattened)
    # plt.matshow(pyem.transition_matrix.transpose())
    # plt.show()

    return trjs, pclims, pcinit, trjs_binned_all0, plot_data



