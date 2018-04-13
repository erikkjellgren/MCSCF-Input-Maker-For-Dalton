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
                print("WARNING: Natural occupation number close to one, might lead to error in choosing active orbitals")
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
                print("WARNING: Natural occupation number close to one, might lead to error in choosing active orbitals")
            if i < 2.0 - threshold and i > 1.5:
                counter_active += 1
            elif i > 2.0 - threshold:
                counter_inactive += 1
        CAS[counter_symmetry] = 2*counter_active
        inactive[counter_symmetry] = counter_inactive
        counter_symmetry += 1
    return CAS, inactive
    