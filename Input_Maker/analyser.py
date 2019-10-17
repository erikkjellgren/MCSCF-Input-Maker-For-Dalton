import numpy as np
import matplotlib.pyplot as plt


def threshold_scan_all(natural_occupations):
    all_orbital = np.array([])
    for key in natural_occupations:
        all_orbital = np.concatenate((all_orbital, natural_occupations[key]))
    all_orbitals = np.sort(all_orbital)
    orbital_counter = 1
    latests_threshold = 0
    weakly_occupied = []
    for i in all_orbitals:
        if i < 1.0:
            weakly_occupied.append(i)
    weakly_occupied = np.array(weakly_occupied)
    weakly_occupied = np.sort(weakly_occupied)[::-1]
    print("Occupation     Correlating    Threshold    Occ-active orb.    Change in threshold")
    for i in all_orbitals:
        if len(weakly_occupied) >= orbital_counter:
            weak_occ = weakly_occupied[orbital_counter-1]
        else:
            weak_occ = 0
        if i != 2.0 and i > 1.0 and latests_threshold != 0:
            print("   {0:.4f}        {1:.4f}        {2:.4f}           {3:2d}                 {4:.4f}".format(i, weak_occ, 2-i, orbital_counter, latests_threshold - (2-i)))
            latests_threshold = 2-i
            orbital_counter += 1
            if orbital_counter == 21:
                break
        elif i != 2.0 and i > 1.0:
            print("   {0:.4f}        {1:.4f}        {2:.4f}           {3:2d}".format(i, weak_occ, 2-i, orbital_counter))
            latests_threshold = 2-i
            orbital_counter += 1
    
    
def threshold_scan_symmetries(natural_occupations):
    all_orbital = np.array([])
    for key in natural_occupations:
        print("Symmetry: ", key)
        orbital_counter = 1
        for i in np.sort(natural_occupations[key]):
            if i != 2.0 and i > 1.0:
                print("Threshold: {0:.4f} gives: {1:2d} occupied-active orbitals".format(2-i, orbital_counter))
                orbital_counter += 1
                if orbital_counter == 21:
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
                
                
def print_relative_natural_occ(relative_natural_occupied, relative_natural_virtuel, natural_occupations, neglect_threshold, show_virtuel):
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
    if show_virtuel == True:
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
                
    
def print_metal_d_orbitals(natural_occupations_unsorted, d_orb, number_occ, number_unocc):
    metal_d_orbitals = []
    for key in natural_occupations_unsorted:
        for i,nat_occ in enumerate(natural_occupations_unsorted[key]):
            metal_d_orbitals.append([nat_occ,d_orb[key][i,0],key,d_orb[key][i,1],d_orb[key][i,2],d_orb[key][i,3],d_orb[key][i,4],d_orb[key][i,5]])
    metal_d_orbitals = np.array(metal_d_orbitals)
    d_orb_idx = np.argsort(metal_d_orbitals[:,0])[::-1]
    occ_idx = d_orb_idx[metal_d_orbitals[d_orb_idx,0]>1.0][::-1]
    unocc_idx = d_orb_idx[metal_d_orbitals[d_orb_idx,0]<1.0]
    occ_d_orbitals = metal_d_orbitals[occ_idx]
    unocc_d_orbitals = metal_d_orbitals[unocc_idx]
    a = 12
    print("Nat. Occ.".ljust(a)+"Orb. #".ljust(a)+"Sym. #".ljust(a)+"dxy/d2-".ljust(a)+"dyz/d1-".ljust(a)+"dzz/d0".ljust(a)+"dxz/d1+".ljust(a)+"dxx-yy/d2+".ljust(a)+"total d")
    for i in range(0, number_occ):
        natocc = occ_d_orbitals[i,0]
        orb = occ_d_orbitals[i,1]
        sym = occ_d_orbitals[i,2]
        dxy = occ_d_orbitals[i,3]
        dyz = occ_d_orbitals[i,4]
        dzz = occ_d_orbitals[i,5]
        dxz = occ_d_orbitals[i,6]
        dxxyy = occ_d_orbitals[i,7]
        total_d = dxy+dyz+dzz+dxz+dxxyy
        print("{0: 0.2f}".format(natocc).ljust(a)+str(int(orb)).ljust(a)+str(int(sym)).ljust(a)+"{0: 0.2f}".format(dxy).ljust(a)+"{0: 0.2f}".format(dyz).ljust(a)+"{0: 0.2f}".format(dzz).ljust(a)+"{0: 0.2f}".format(dxz).ljust(a)+"{0: 0.2f}".format(dxxyy).ljust(a)+"{0: 0.2f}".format(total_d))
    print("======================================================================================================")
    for i in range(0, number_unocc):
        natocc = unocc_d_orbitals[i,0]
        orb = unocc_d_orbitals[i,1]
        sym = unocc_d_orbitals[i,2]
        dxy = unocc_d_orbitals[i,3]
        dyz = unocc_d_orbitals[i,4]
        dzz = unocc_d_orbitals[i,5]
        dxz = unocc_d_orbitals[i,6]
        dxxyy = unocc_d_orbitals[i,7]
        total_d = dxy+dyz+dzz+dxz+dxxyy
        print("{0: 0.2f}".format(natocc).ljust(a)+str(int(orb)).ljust(a)+str(int(sym)).ljust(a)+"{0: 0.2f}".format(dxy).ljust(a)+"{0: 0.2f}".format(dyz).ljust(a)+"{0: 0.2f}".format(dzz).ljust(a)+"{0: 0.2f}".format(dxz).ljust(a)+"{0: 0.2f}".format(dxxyy).ljust(a)+"{0: 0.2f}".format(total_d))


def print_metal_d_orbitals_hf(hf_orb_energies, d_orb, number_closed_shells, number_occ, number_unocc):
    metal_d_orbitals_occ = []
    metal_d_orbitals_unocc = []
    counter = 0
    for key in hf_orb_energies:
        for i in range(0, int(number_closed_shells[counter])):
            metal_d_orbitals_occ.append([hf_orb_energies[key][i],d_orb[key][i,0],key,d_orb[key][i,1],d_orb[key][i,2],d_orb[key][i,3],d_orb[key][i,4],d_orb[key][i,5]])
        counter += 1
    counter = 0
    for key in hf_orb_energies:
        for i in range(int(number_closed_shells[counter]), len(hf_orb_energies[key])):
            metal_d_orbitals_unocc.append([hf_orb_energies[key][i],d_orb[key][i,0],key,d_orb[key][i,1],d_orb[key][i,2],d_orb[key][i,3],d_orb[key][i,4],d_orb[key][i,5]])
        counter += 1
    metal_d_orbitals_occ = np.array(metal_d_orbitals_occ)
    metal_d_orbitals_unocc = np.array(metal_d_orbitals_unocc)
    occ_idx = np.argsort(metal_d_orbitals_occ[:,0])[::-1]
    unocc_idx = np.argsort(metal_d_orbitals_unocc[:,0])
    occ_d_orbitals = metal_d_orbitals_occ[occ_idx]
    unocc_d_orbitals = metal_d_orbitals_unocc[unocc_idx]
    a = 12
    print("Orb. en.".ljust(a)+"Orb. #".ljust(a)+"Sym. #".ljust(a)+"dxy/d2-".ljust(a)+"dyz/d1-".ljust(a)+"dzz/d0".ljust(a)+"dxz/d1+".ljust(a)+"dxx-yy/d2+".ljust(a)+"total d")
    for i in range(0, number_occ):
        orb_en = occ_d_orbitals[i,0]
        orb = occ_d_orbitals[i,1]
        sym = occ_d_orbitals[i,2]
        dxy = occ_d_orbitals[i,3]
        dyz = occ_d_orbitals[i,4]
        dzz = occ_d_orbitals[i,5]
        dxz = occ_d_orbitals[i,6]
        dxxyy = occ_d_orbitals[i,7]
        total_d = dxy+dyz+dzz+dxz+dxxyy
        print("{0: 0.2E}".format(orb_en).ljust(a)+str(int(orb)).ljust(a)+str(int(sym)).ljust(a)+"{0: 0.2f}".format(dxy).ljust(a)+"{0: 0.2f}".format(dyz).ljust(a)+"{0: 0.2f}".format(dzz).ljust(a)+"{0: 0.2f}".format(dxz).ljust(a)+"{0: 0.2f}".format(dxxyy).ljust(a)+"{0: 0.2f}".format(total_d))
    print("======================================================================================================")
    for i in range(0, number_unocc):
        orb_en = unocc_d_orbitals[i,0]
        orb = unocc_d_orbitals[i,1]
        sym = unocc_d_orbitals[i,2]
        dxy = unocc_d_orbitals[i,3]
        dyz = unocc_d_orbitals[i,4]
        dzz = unocc_d_orbitals[i,5]
        dxz = unocc_d_orbitals[i,6]
        dxxyy = unocc_d_orbitals[i,7]
        total_d = dxy+dyz+dzz+dxz+dxxyy
        print("{0: 0.2E}".format(orb_en).ljust(a)+str(int(orb)).ljust(a)+str(int(sym)).ljust(a)+"{0: 0.2f}".format(dxy).ljust(a)+"{0: 0.2f}".format(dyz).ljust(a)+"{0: 0.2f}".format(dzz).ljust(a)+"{0: 0.2f}".format(dxz).ljust(a)+"{0: 0.2f}".format(dxxyy).ljust(a)+"{0: 0.2f}".format(total_d))









