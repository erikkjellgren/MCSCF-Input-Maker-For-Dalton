import numpy as np


def Pick_RAS_active_threshold(threshold, Natural_Occupations):
    RAS1 = np.zeros(len(Natural_Occupations), dtype=int)
    RAS3 = np.zeros(len(Natural_Occupations), dtype=int)
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    counter_symmetry = 0
    for key in Natural_Occupations:
        counter_active = 0
        counter_inactive = 0
        for i in Natural_Occupations[key]:
            if i < 1.5 and i > 0.5:
                print("WARNING: Natural occupation number close to one, might lead to error in choosing active orbitals.")
            if i < 2.0 - threshold and i > 1.4:
                counter_active += 1
            elif i > 2.0 - threshold:
                counter_inactive += 1
        RAS1[counter_symmetry] = counter_active
        RAS3[counter_symmetry] = counter_active
        inactive[counter_symmetry] = counter_inactive
        counter_symmetry += 1
    return RAS1, RAS3, inactive
    
    
def Pick_CAS_active_threshold(threshold, Natural_Occupations):
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    counter_symmetry = 0
    for key in Natural_Occupations:
        counter_active = 0
        counter_inactive = 0
        for i in Natural_Occupations[key]:
            if i < 1.5 and i > 0.5:
                print("WARNING: Natural occupation number close to one, might lead to error in choosing active orbitals.")
            if i < 2.0 - threshold and i > 1.5:
                counter_active += 1
            elif i > 2.0 - threshold:
                counter_inactive += 1
        CAS[counter_symmetry] = 2*counter_active
        inactive[counter_symmetry] = counter_inactive
        counter_symmetry += 1
    return CAS, inactive
    

def Pick_CAS_number_occupied(number_occ, Natural_Occupations, allow_more_virt=False):
    """
    Function to pick a CAS based on number of wanted occupied active orbitals.
    The number of virtuel orbitals will be the same as the number of occupied.
    The occupied and virtuel is picked in the order, of being most different,
      from 2 or 0 respectively.
    If number_occ is larger than what is possible, extra vitual orbitals can still 
      be added if allow_more_virt=True is set
    """
    number_electrons = 0.5
    for key in Natural_Occupations:
        number_electrons += np.sum(Natural_Occupations[key])
    number_electrons = int(number_electrons)
    # Incase allow_more_virt=True, number_virt can be larger than number_occ
    number_virt = number_occ 
    # Cannot choose more occupied orbitals than allowed by number of electrons
    number_occ = np.min([number_occ, number_electrons//2])
    if allow_more_virt == False:
        # Ensure number of virt == number of occ
        number_virt = number_occ
        
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    # array to hold symmetry and occupation numbers of picked orbitals
    picked_occ = np.zeros((number_occ,2))
    picked_occ[:,:] = 3 # Set to three, all other occupations will be smaller
    picked_virt = np.zeros((number_virt,2))
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    for key in Natural_Occupations:
        for i in Natural_Occupations[key]:
            if i < 1.1 and i > 0.9:
                print("WARNING: Cannot determine if orbital is occupied or virtuel, might lead to error in choosing orbitals.")
                if i > 1.0:
                    print("Occupation: "+str(i)+" chosen as occupied.")
                else:
                    print("Occupation: "+str(i)+" chosen as virtuel.")
            if i > 1.0:
                inactive[key-1] += 1 # Will subtract picked orbitals later
                j = np.argmax(picked_occ[:,1])
                if i < picked_occ[j,1]:  
                    picked_occ[j,0] = key
                    picked_occ[j,1] = i
            else:
                j = np.argmin(picked_virt[:,1])
                if i > picked_virt[j,1]:  
                    picked_virt[j,0] = key
                    picked_virt[j,1] = i
    for i in range(len(picked_occ[:,0])):
        idx = int(picked_occ[i,0]) - 1
        if picked_occ[i,1] != 2:
            CAS[idx] += 1
            inactive[idx] -= 1
    for i in picked_virt[:,0]:
        CAS[int(i)-1] += 1
    return CAS, inactive


def Pick_RASCI_number_occupied(number_occ, Natural_Occupations, max_virtuel, excitation=[0,0]):
    """
    Function to pick a RAS-CI active based on number number of 
      occupied orbitals wanted. 
      
    excitation = [from symmetry, to symmetry]
    """
    number_electrons = 0.5
    for key in Natural_Occupations:
        number_electrons += np.sum(Natural_Occupations[key])
    number_electrons = int(number_electrons)
    # Cannot choose more occupied orbitals than allowed by number of electrons
    number_occ = np.min([number_occ, number_electrons//2])
    RAS1 = np.zeros(len(Natural_Occupations), dtype=int)
    RAS2 = np.zeros(len(Natural_Occupations), dtype=int)
    RAS2_electrons = 0
    RAS3 = np.zeros(len(Natural_Occupations), dtype=int)
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    # array to hold symmetry and occupation numbers of picked orbitals
    picked_occ = np.zeros((number_occ,2))
    picked_occ[:,:] = 3 # Set to three, all other occupations will be smaller
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    for key in Natural_Occupations:
        for i in Natural_Occupations[key]:
            if i < 1.1 and i > 0.9:
                print("WARNING: Cannot determine if orbital is occupied or virtuel, might lead to error in choosing orbitals.")
                if i > 1.0:
                    print("Occupation: "+str(i)+" chosen as occupied.")
                else:
                    print("Occupation: "+str(i)+" chosen as virtuel.")
            if i > 1.0:
                inactive[key-1] += 1 # Will subtract picked orbitals later
                j = np.argmax(picked_occ[:,1])
                if i < picked_occ[j,1]:
                    picked_occ[j,0] = key
                    picked_occ[j,1] = i
    for i in range(len(picked_occ[:,0])):
        idx = int(picked_occ[i,0]) - 1
        if picked_occ[i,1] != 2:
            RAS1[idx] += 1
            inactive[idx] -= 1
            
    for key in Natural_Occupations:
        nat_occ_temp = np.zeros((len(Natural_Occupations[key]),2))
        nat_occ_temp[:,0] = Natural_Occupations[key]
        nat_occ_temp[:,1] = key
        if key == 1:
            nat_occ_all = nat_occ_temp
        else:
            nat_occ_all = np.vstack((nat_occ_all, nat_occ_temp))
    idx = np.argsort(nat_occ_all[:,0])[::-1]
    nat_occ_all[:,0] = nat_occ_all[idx,0]
    nat_occ_all[:,1] = nat_occ_all[idx,1]
    for i in range(0, len(nat_occ_all)):
        if np.sum(RAS3) == max_virtuel:
            break
        nat_occ = nat_occ_all[i,0]
        symmetry = nat_occ_all[i,1]
        if nat_occ < 1:
            RAS3[int(symmetry)-1] += 1
    if excitation != [0,0]:
        if RAS1[excitation[0]-1] != 0 and RAS3[excitation[1]-1] != 0:
            RAS1[excitation[0]-1] -= 1
            RAS2[excitation[0]-1] += 1
            RAS2_electrons += 2
            RAS3[excitation[1]-1] -= 1
            RAS2[excitation[1]-1] += 1
    return RAS1, RAS2, RAS3, inactive, RAS2_electrons


def Pick_CAS_threshold_electron_retrieval(occupied_threshold, electron_retrieval, Natural_Occupations, number_of_symmetries, print_retrieved_electron):
    """
    Heuristic to pick CAS.
    
    Input : occupied_threshold, denotes a threshold for which occupied
            orbitals should be included in CAS
            If occupied_threshold = 1.98, then all occupied orbitals with
            a natural occupation less than 1.98 will be included (and an
            occupation above 1.0).
            electron_retrieval, denotes how many virtuel orbitals should be
            included. If electron_retrieval=0.8, then virtuel orbitals will
            be included untill 80% of the electron missing from the picked
            occupied is retrieved.
            
    The choosing of virtuel orbitals is symmetry independent.
    
    Assume that the keys are the numbers from 1 and upwards
    """
    CAS = np.zeros(number_of_symmetries, dtype=int)
    inactive = np.zeros(number_of_symmetries, dtype=int)
    missing_electron = 0 # Counts how different the occupied orbitals are from zero
    for key in range(1, number_of_symmetries+1):
        counter_active = 0
        counter_inactive = 0
        if key in Natural_Occupations:
            for i in Natural_Occupations[key]:
                if i < 1.1 and i > 0.9:
                    print("WARNING: Natural occupation number close to one, might lead to error in choosing active orbitals.")
                if i >= 1.0 and i <= occupied_threshold:
                    counter_active += 1
                    missing_electron += 2 - i
                elif i > occupied_threshold:
                    counter_inactive += 1
        CAS[key-1] = counter_active
        inactive[key-1] = counter_inactive
    if missing_electron == 0.0:
        # If no electron is missing, then no occupied orbitals
        # have been added, therefore no virtuel should be added
        return CAS, inactive
    # Pick virtuel orbitals from here
    all_virtuel = []
    all_virtuel_sym = []
    for key in Natural_Occupations:
        # Make a list of all virtuel and their symmetry
        for occupation in Natural_Occupations[key]:
            if occupation < 1.0:
                all_virtuel.append(occupation)
                all_virtuel_sym.append(key)
    electron_retrieved = 0.0
    for i in range(0, 30):
        # Pick the largest occupation and then zero it
        index = np.argmax(all_virtuel)
        CAS[all_virtuel_sym[index]-1] += 1
        electron_retrieved += all_virtuel[index]
        electron_retrieved_percentage = 1 - (missing_electron - electron_retrieved)/missing_electron
        all_virtuel[index] = 0.0
        if electron_retrieved_percentage > electron_retrieval:
            break
    if print_retrieved_electron == True:
        print("Percentage electron retrieved: {0:0.2f}".format(electron_retrieved_percentage))
    return CAS, inactive




