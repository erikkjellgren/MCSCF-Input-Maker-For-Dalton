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
    

def Pick_CAS_number_occupied(number_occ, Natural_Occupations):
    """
    Function to pick a CAS based on number of wanted occupied active orbitals.
    The number of virtuel orbitals will be the same as the number of occupied.
    The occupied and virtuel is picked in the order, of being most different,
      from 2 or 0 respectively.
    """
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    # array to hold symmetry and occupation numbers of picked orbitals
    picked_occ = np.zeros((number_occ,2))
    picked_occ[:,:] = 2 # Set to two, all other occupations will be smaller
    picked_virt = np.zeros((number_occ,2))
    inactive = np.zeros(len(Natural_Occupations), dtype=int)
    CAS = np.zeros(len(Natural_Occupations), dtype=int)
    for key in Natural_Occupations:
        for i in Natural_Occupations[key]:
            if i < 1.1 and i > 0.9:
                print("WARNING: Cannot determine if orbital is occupied or virtuel, might lead to error in choosing orbitals.")
            if i > 1.1:
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
    for i in picked_occ[:,0]:
        CAS[int(i)-1] += 1
        inactive[int(i)-1] -= 1
    for i in picked_virt[:,0]:
        CAS[int(i)-1] += 1
    return CAS, inactive
    







