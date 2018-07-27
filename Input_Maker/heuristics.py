import numpy as np
from scipy.special import binom


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


def Pick_RASCI_number_occupied(number_occ, Natural_Occupations, approx_determinant_limt=7*10**6, excitation=[0,0]):
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
        for i in Natural_Occupations[key]:
            if i < 1.0:
                if np.sum(RAS3)*2+2 - np.sum(RAS1)*2 > 1:
                    if binom(np.sum(RAS1)*2,2)*binom(np.sum(RAS3)*2+2 - np.sum(RAS1)*2,2) < approx_determinant_limt:
                        RAS3[key-1] += 1
                    else:
                        break
                else:
                    RAS3[key-1] += 1
    if excitation != [0,0]:
        if RAS1[excitation[0]-1] != 0 and RAS3[excitation[1]-1] != 0:
            RAS1[excitation[0]-1] -= 1
            RAS2[excitation[0]-1] += 1
            RAS2_electrons += 2
            RAS3[excitation[1]-1] -= 1
            RAS2[excitation[1]-1] += 1

    return RAS1, RAS2, RAS3, inactive, RAS2_electrons







