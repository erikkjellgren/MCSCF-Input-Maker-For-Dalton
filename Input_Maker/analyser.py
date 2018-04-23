import numpy as np
import matplotlib.pyplot as plt


def threshold_scan_all(natural_occupations):
    all_orbital = np.array([])
    for key in natural_occupations:
        all_orbital = np.concatenate((all_orbital, natural_occupations[key]))
    all_orbitals = np.sort(all_orbital)
    orbital_counter = 1
    latests_threshold = 0
    print("Threshold    Occ-active orb.    Change in threshold")
    for i in all_orbitals:
        if i != 2.0 and i > 1.4 and latests_threshold != 0:
            print("  {0:.3f}            {1:2d}                 {2:.3f}".format(2-i, orbital_counter, latests_threshold - (2-i)))
            latests_threshold = 2-i
            orbital_counter += 1
            if orbital_counter == 17:
                break
        elif i != 2.0 and i > 1.4:
            print("  {0:.3f}            {1:2d}".format(2-i, orbital_counter))
            latests_threshold = 2-i
            orbital_counter += 1
    
    
def threshold_scan_symmetries(natural_occupations):
    all_orbital = np.array([])
    for key in natural_occupations:
        print("Symmetry: ", key)
        orbital_counter = 1
        for i in np.sort(natural_occupations[key]):
            if i != 2.0 and i > 1.4:
                print("Threshold: {0:.4f} gives: {1:2d} occupied-active orbitals".format(2-i, orbital_counter))
                orbital_counter += 1
                if orbital_counter == 17:
                    break


def print_natural_occ(natural_occupations, threshold, neglect_threshold):
    for key in natural_occupations:
        if key != 1:
            print("")
            print("")
        print("Symmetry: ", key)
        orbitals = np.resize(natural_occupations[key][(natural_occupations[key]>neglect_threshold) & (natural_occupations[key]<2-neglect_threshold)], ((len(natural_occupations[key][(natural_occupations[key]>neglect_threshold) & (natural_occupations[key]<2-neglect_threshold)])+1)//5,5))
        latest = 2
        for i in range(len(orbitals)):
            for j in range(len(orbitals[0])):
                if orbitals[i,j] > latest:
                    orbitals[i,j] = np.nan
                else:
                    latest = orbitals[i,j]

        counter = 0
        for i in range(len(orbitals)):
            print("")
            for j in range(len(orbitals[0])):
                if 2 - orbitals[i,j] > threshold and orbitals[i,j] > 0.5:
                    counter += 1
                    print('\033[4m' + "{0:.4f}".format(orbitals[i,j])+ '\033[0m',end="")
                    print("   ", end="")
                elif orbitals[i,j] < 0.5 and counter != 0:
                    print('\033[4m' + "{0:.4f}".format(orbitals[i,j])+ '\033[0m',end="")
                    print("   ", end="")
                    counter -= 1
                else:
                    if np.isnan(orbitals[i,j]):
                        print("         ",end="")
                    else:
                        print("{0:.4f}   ".format(orbitals[i,j]),end="")