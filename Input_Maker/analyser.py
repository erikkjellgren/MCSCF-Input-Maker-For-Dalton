import numpy as np
import matplotlib.pyplot as plt


def threshold_scan_all(natural_occupations):
    all_orbital = np.array([])
    for key in natural_occupations:
        all_orbital = np.concatenate((all_orbital, natural_occupations[key]))
    all_orbitals = np.sort(all_orbital)
    orbital_counter = 1
    for i in all_orbitals:
        if i != 2.0 and i > 1.4:
            print("Threshold: {0:.4f} gives: {1:2d} occupied-active orbitals".format(2-i, orbital_counter))
            orbital_counter += 1
            if orbital_counter == 17:
                break
    
    
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
        A1 = np.resize(natural_occupations[key][(natural_occupations[key]>neglect_threshold) & (natural_occupations[key]<2-neglect_threshold)], ((len(natural_occupations[key][(natural_occupations[key]>neglect_threshold) & (natural_occupations[key]<2-neglect_threshold)])+1)//5,5))
        counter = 0
        for i in range(len(A1)):
            for j in range(len(A1[0])):
                counter += 1
                if counter > len(natural_occupations[key]):
                    A1[i,j] = np.nan
        counter = 0
        for i in range(len(A1)):
            print("")
            for j in range(len(A1[0])):
                if 2 - A1[i,j] > threshold and A1[i,j] > 0.5:
                    counter += 1
                    print('\033[4m' + "{0:.4f}".format(A1[i,j])+ '\033[0m',end="")
                    print("   ", end="")
                elif A1[i,j] < 0.5 and counter != 0:
                    print('\033[4m' + "{0:.4f}".format(A1[i,j])+ '\033[0m',end="")
                    print("   ", end="")
                    counter -= 1
                else:
                    if np.isnan(A1[i,j]):
                        print("         ",end="")
                    else:
                        print("{0:.4f}   ".format(A1[i,j]),end="")