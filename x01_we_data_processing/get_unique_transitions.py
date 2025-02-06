#-----------------------------------------------------------------------------------------------------------------------
#Jonathan Borowsky
#5/7/2024
#Grabe lab
#                   Mature molecular dynamics data processing methods from jupyter notebooks
#copied from westpa_msm_functions.py to avoid pyemma import issues; it appears not to be compatible with the lastest version of python
#-----------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------
#                                        MACROSTATE TRANSITION ANALYSIS
#          determine the number of independent transitions between macrostates in a westpa simulation
#                                    from westpa-barrier-crossing-v2.ipynb
#-----------------------------------------------------------------------------------------------------------------------

import numpy as np
import h5py

# ------------------------
# Parameters:
# h5path: path to the input west.h5 file
# pc_2_macrostate: function which takes a progress coordinate as its sole required argument
#                  and returns an integer in the range [0,n_macrostates) indicating which
#                  macrostate the walker is in, with -1 indicating no macrostate
# n_macrostates: number of macrostates in the system
# minround: first round of the h5 file to use. values other than 0 not currently supported
# maxround: last round of the h5 file to use. -1 means to use all rounds
# allcomplete: False means not to use the last round in the h5 file because that file
#              contains data for only some walkers from an incomplete westpa round
# ------------------------
# Returns:
# transitions: westpa rounds and walker numbers of walkers which have reached each macrostate from each other macrostate
#              If a walker leaves one macrostate and splits, and multiple of its offspring reach a second macrostate,
#              only the first of these offspring to reach the second state is recorded
#              (treating lower walker number as earlier within each westpa round)

def get_unique_transitions(h5path, pc_2_macrostate, n_macrostates=2, minround=0, maxround=-1, allcomplete=False, printfreq=10):
    # each entry describes:
    # which macrostate has has the trajectory leading to each walker most recently been in, 
    # which ancestor was most recently in that state,
    # whether the walker is currently in either macrostate,
    # the most recent ancestor to have reached the current macrostate from elsewhere
    ensembles = []

    # list of all walkers with at least one series of descendants
    # which reached each state from each other state without returning (or remaining any longer in their initial state)
    transition_start_states = [[[] for i in range(n_macrostates)] for i in range(n_macrostates)]

    pcs_by_tss = {}
    
    # list of walkers which have transitioned from state i to j without passing through any other macrostate,
    # each from a unique parent
    # for counting transitions and knowing which walkers to start reconstructing transition paths from
    transitions = [[[] for i in range(n_macrostates)] for i in range(n_macrostates)]

    #walkers_open = 0
    walkers_by_founder = {}
    independent_transition_timings = []
    
    # load h5 file
    with h5py.File(h5path, 'r') as f:

        # determine total number of westpa rounds and set number to use accordingly
        iterations = [it for it in f["iterations"]]
        maxiter = len(iterations)
        if maxround == -1:
            if allcomplete:
                maxround = maxiter
            else:  # if the last round is not complete, omit it as data have yet to be written
                maxround = maxiter - 1

        for i in range(minround, maxround):
            if i % printfreq == 0:
                print(f"round {i}")

            # this extracts only the iteration name, not the iteration data
            iter_name = iterations[i]
            # using the iteration name to extract the data
            iter_data = f["iterations"][iter_name]

            if i == 0:
                # progress coordinate of the starting structure (0th component for multidimensional PCs)
                init_pc0 = [i for i in iter_data["ibstates"]["bstate_pcoord"]][0]

                init_macrostate = pc_2_macrostate(init_pc0)

                if init_macrostate == -1:
                    print("error; system begins between macrostates")
                    return 0

                ensembles.append([[init_macrostate, -1, 0, 0, -1, 0]])

            ensemble = []

            for j, (pc_, parent_) in enumerate(zip(iter_data["pcoord"], iter_data["seg_index"])):

                pc = pc_[0]
                parent = parent_[1]

                if i == 0 and parent < 0:
                    parent = -1

                # determine which if any macrostate the walker is currently in
                current_state = pc_2_macrostate(pc)

                # if the walker is in a macrostate
                if current_state != -1:

                    #for branch inspection, not used to produce output
                    # if (ensembles[-1][parent][0] != current_state and current_state != 0) \
                    #         and ensembles[-1][parent][1:3] in \
                    #         transition_start_states[ensembles[-1][parent][0]][current_state]:
                    #     print(f"walker {i},{j} reached state {current_state} but another relative descended from the same ancestor last in state  {ensembles[-1][parent][0]} has already reached {current_state}")
                    
                    # if the parent of this walker left a macrostate
                    # and its other descendants have not already reached this one
                    # [without passing through a third macrostate], record a transition
                    if (ensembles[-1][parent][0] != current_state or ensembles[-1][parent][3] == 1) \
                            and ensembles[-1][parent][1:3] not in \
                            transition_start_states[ensembles[-1][parent][0]][current_state]:
                        
                        transitions[ensembles[-1][parent][0]][current_state].append([i, j])
                        transition_start_states[ensembles[-1][parent][0]][current_state].append(
                            ensembles[-1][parent][1:3])


                        if current_state != 0 and ensembles[-1][parent][0] != current_state:
                            #print(f"{i},{j}, pc={pc}")
                            
                            #PC_TSS
                            pcs_by_tss[repr(ensembles[-1][parent][1:3])] = [[[i,j]], [[i,j,pc[0]]]]
                        
                    
                        #record this trajectory as the founder of a new subpopulation of walkers in this macrostate
                        # only if it originates from a different macrostate; 
                        # descendants of recrossers are still part of their original population
                        if ensembles[-1][parent][0] != current_state:
                            founder = [i,j]
                            if current_state != 0:
                                independent_transition_timings.append(i)
                        else:
                            founder = ensembles[-1][parent][4:6]

                    #this covers the case where a walker departs from state A, splits in between macrostates, and has multiple descendants reach state B. Only the first such descendant is counted as an independent transition (as only one such transition has occurred), but all of them are counted as independent founders of subpopulations in state B. Part of the reason for this is that just taking the first descendant to reach B is arbitrary and will not give you an accurate picture of the set of substates being sampled in B.
                    elif ensembles[-1][parent][0] != current_state:
                        founder = [i,j]
                        
                        #PC_TSS
                        if current_state != 0:
                            pcs_by_tss[repr(ensembles[-1][parent][1:3])][0].append([i,j])
                            pcs_by_tss[repr(ensembles[-1][parent][1:3])][1].append([i,j,pc[0]])
                        
                    else:
                        #if this walker is not newly arrived from a different macrostate, the most recent ancestor to reach the current macrostate 
                        # from elsewhere is the same as for this walker's parent
                        founder = ensembles[-1][parent][4:6]
                        
                        #PC_TSS
                        if current_state != 0:
                            for k in pcs_by_tss.keys():
                                if founder in pcs_by_tss[k][0]:
                                    pcs_by_tss[k][1].append([i,j,pc[0]])
                                    
                            # pcs_by_tss[ensembles[-1][parent][1:3]][0].append([i,j])
                            # pcs_by_tss[ensembles[-1][parent][1:3]][1].append([i,j,pc])
                        
                    # record the round and walker at which the trajectory
                    # was most recently in its present macrostate
                    ensemble.append([current_state, i, j, 0] + founder)
                    
                    if current_state != 0:
                        wr = repr(founder + [current_state])
                        if wr not in walkers_by_founder.keys():
                            walkers_by_founder[wr] = []
                        
                        walkers_by_founder[wr].append([i,j])
                    
                else:
                    # if the walker is not in any macrostate,
                    # its most recent ancestor in either macrostate is the same as that of its parent
                    # (and in some cases may be its parent)
                    ensemble.append(ensembles[-1][parent][0:3] + [1] + ensembles[-1][parent][4:6])

            ensembles.append(ensemble)

    #print(walkers_open)
    return transitions, walkers_by_founder, independent_transition_timings, pcs_by_tss
