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
    """
    Function to print natural occupations in a formatted way.
    
    Input : natural_occupations, dictionary of ordered natural occupations.
          : threshold, threshold for which occupation numbers will be underlined.
               Threshold is set for how close to 2, based on this the same number of
               occupied and unoccupied will be underlined in each symmetry.
          : neglect_threshold, thershold for which occupation numbers will not be printed
    """
    for key in natural_occupations:
        if key != 1:
            print("")
            print("")
        print("Symmetry: ", key)
        counter = 0
        orbital_counter = 0
        for i in range(len(natural_occupations[key])):
            if natural_occupations[key][i] < neglect_threshold or natural_occupations[key][i] > 2-neglect_threshold:
                continue
            if orbital_counter%5 == 0:
                print("")
            orbital_counter += 1
            if 2 - natural_occupations[key][i] > threshold and natural_occupations[key][i] > 0.5:
                counter += 1
                print('\033[4m' + "{0:.4f}".format(natural_occupations[key][i])+ '\033[0m',end="")
                print("   ", end="")
            elif natural_occupations[key][i] < 0.5 and counter != 0:
                print('\033[4m' + "{0:.4f}".format(natural_occupations[key][i])+ '\033[0m',end="")
                print("   ", end="")
                counter -= 1
            else:
                print("{0:.4f}   ".format(natural_occupations[key][i]),end="")
                
                
def print_relative_natural_occ(relative_natural_occupied, relative_natural_virtuel, natural_occupations, neglect_threshold):
    """
    Function to print relative natural occupations.
    
    natural_occupations is used to make sure neglect_threshold will make the same
      orbitals show for normal and relative print.
    """
    for key in relative_natural_occupied:
        if key != 1:
            print("")
            print("")
        print("Symmetry (occupied): ", key)
        orbital_counter = 0
        for i in range(len(relative_natural_occupied[key])):
            if natural_occupations[key][i] > 2-neglect_threshold:
                continue
            if orbital_counter%5 == 0:
                print("")
            orbital_counter += 1
            print("{0:.4f}   ".format(relative_natural_occupied[key][i]),end="")
    for key in relative_natural_virtuel:
        print("")
        print("")
        print("Symmetry (virtuel): ", key)
        orbital_counter = 0
        for i in range(len(relative_natural_virtuel[key])):
            # len(relative_natural_occupied[key]) to find virtuel orbitals in the full list
            if natural_occupations[key][i+len(relative_natural_occupied[key])] < neglect_threshold:
                continue
            if orbital_counter%5 == 0:
                print("")
            orbital_counter += 1
            print("{0:.4f}   ".format(relative_natural_virtuel[key][i]),end="")
                
    











