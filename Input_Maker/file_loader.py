import numpy as np


def wavefunction_type_output(load_file):
    found_check = 0
    for i in load_file:
        if i[0:23] == "@    Wave function type":
            data = i[24:].split()
            for j in data:
                if j != '' and j != "---":
                    wavefunction = j.lower()
                    found_check = 1
                    break
    if found_check == 0:
        print("Warning: Wave function type could not be determined. Set to \"Unknown\".")
        print("This might cause problems if it is a MP2 wave function.")
        wavefunction = "Unknown"
    return wavefunction


def total_nuclei_charge(load_file):
    for i in load_file:
        if i[0:8] == "  total:":
            data = i[9:].split()
            number_counter = 0
            for j in data:
                if j != '':
                    number_counter += 1
                if number_counter == 2:
                    charge = int(float(j))
                    break
    return charge


def orbital_symmetries(load_file):
    orbital_symmetries = np.zeros(20)
    for i in load_file:
        if i[0:38] == "  Number of orbitals in each symmetry:":
            data = i[40:].split(" ")
            counter = 0
            for j in data:
                if j != '':
                    orbital_symmetries[counter] = int(j)
                    counter += 1
            break
        else:
            counter = 1
    return orbital_symmetries[0:counter]
    
    
def closed_shell_number(load_file):
    closed_shells = np.zeros(9)
    for i in load_file:
        if i[0:26] == "    Closed shell orbitals:":
            data = i[27:].split(" ")
            counter = 0
            for j in data:
                if j != '':
                    closed_shells[counter] = int(j)
                    counter += 1
            break
    return closed_shells[0:counter]
    
    
def closed_shell_number_hf(load_file):
    closed_shells = np.zeros(9)
    for line in load_file:
        if "Orbital occupations :" in line:
            data = line.split()
            counter = 0
            for j in range(3, len(data)):
                closed_shells[counter] = int(data[j])
                counter += 1
            break
    return closed_shells[0:counter]
    

def electronsMP2(load_file):
    found_electrons = 0 # Check if number electrons is found
    electrons = 0
    for i in load_file:
        if i[0:26] == "    Number of electrons  :":
            data = i[27:].split(' ')
            for j in data:
                if j != '':
                    electrons = int(j)
                    found_electrons = 1
                    break
            if found_electrons == 1:
                break
    return electrons
    
    
def electrons(load_file):
    found_electrons = 0 # Check if number electrons is found
    electrons = 0
    for i in load_file:
        if i[0:37] == "@    Number of closed shell electrons":
            data = i[38:].split(' ')
            for j in data:
                if j != '':
                    electrons += int(j)
                    break
        if i[0:41] == "@    Number of electrons in active shells":
            data = i[42:].split(' ')
            for j in data:
                if j != '':
                    electrons += int(j)
                    found_electrons = 1
                    break
            if found_electrons == 1:
                break
    return electrons
    
    
def HF_orb_energies(load_file, number_symmetries):
    HF_energies = {}
    get_energies_check = 0 # Check when HF energies are found
    energies_string = '' # Accumalte HF energies
    energies_array = np.zeros(1000) # Array to temp store energies
    newline_check = 0 # Check newlines for stop
    current_symmetry = 1
    for i in load_file:
        if i[0:42] == " Hartree-Fock orbital energies, symmetry "+str(current_symmetry):
            get_energies_check = 1
        elif get_energies_check == 1:
            if i == "\n":
                newline_check += 1
            energies_string = energies_string + i
        if newline_check == 2:
            get_energies_check = 0 # Reset
            newline_check = 0 # Reset
            
            data = energies_string.split(" ")
            counter = 0
            for j in data:
                if j != '' and j != "\n":
                    energies_array[counter] = float(j)
                    counter += 1
            HF_energies[current_symmetry] = np.copy(energies_array[0:counter])
            energies_string = '' # Reset
            current_symmetry += 1 # Go to next symmetry
        if current_symmetry  == number_symmetries+1:
            break         
    return HF_energies
    
    
def HF_orb_energies_hf_wf(load_file, number_symmetries):
    HF_energies = {}
    get_energies_check = False # Check when HF energies are found
    current_symmetry = 0
    for line in load_file:
        if " E(LUMO)" in line:
            HF_energies[current_symmetry] = np.array(energies)
            break   
        if "Hartree-Fock orbital energies" in line:
            get_energies_check = True
        elif len(line) > 10 and get_energies_check:
            orb_energies = line.split()
            if len(orb_energies) == 7:
                if current_symmetry != 0:
                    HF_energies[current_symmetry] = np.array(energies)
                energies = []
                current_symmetry += 1 # Go to next symmetry
                for i in range(2, 7):
                    energies.append(float(orb_energies[i]))
            else:
                for en in orb_energies:
                    energies.append(float(en))
    return HF_energies
    

def Natural_Occupations_MP2(load_file, number_symmetries):
    Natural_Occupations = {}
    get_occupations_check = 0 # Check when HF energies are found
    occupations_string = '' # Accumalte HF energies
    occupations_array = np.zeros(1000) # Array to temp store energies
    newline_check = 0 # Check newlines for stop
    current_symmetry = 1
    for i in load_file:
        if i[0:47] == " Natural orbital occupation numbers, symmetry "+str(current_symmetry):
            get_occupations_check = 1
        elif get_occupations_check == 1:
            if i == "\n":
                newline_check += 1
            if newline_check == 1:
                occupations_string = occupations_string + i
        if newline_check == 2:
            get_occupations_check = 0 # Reset
            newline_check = 0 # Reset
            
            data = occupations_string.split(" ")
            counter = 0
            for j in data:
                if j != '' and j != "\n":
                    occupations_array[counter] = float(j)
                    counter += 1
            Natural_Occupations[current_symmetry] = np.copy(occupations_array[0:counter])
            occupations_string = '' # Reset
            current_symmetry += 1 # Go to next symmetry
        if current_symmetry  == number_symmetries+1:
            break
            
    return Natural_Occupations
    
    
def Natural_Occupations(load_file, number_symmetries):
    Natural_Occupations = {}
    get_occupations_check = 0 # Check when HF energies are found
    occupations_string = '' # Accumalte HF energies
    occupations_array = np.zeros(1000) # Array to temp store energies
    newline_check = 0 # Check newlines for stop
    current_symmetry = 1
    for i in load_file:
        if i[0:11] == " Symmetry "+str(current_symmetry):
            get_occupations_check = 1
            if "No occupied orbitals" in i:
                get_occupations_check = 0
                current_symmetry += 1
        elif get_occupations_check == 1:
            if i == "\n":
                newline_check += 1
            occupations_string = occupations_string + i
        if newline_check == 2:
            get_occupations_check = 0 # Reset
            newline_check = 0 # Reset
            
            data = occupations_string.split(" ")
            counter = 0
            for j in data:
                if j != '' and j != "\n":
                    occupations_array[counter] = float(j)
                    counter += 1
            Natural_Occupations[current_symmetry] = np.copy(occupations_array[0:counter])
            occupations_string = '' # Reset
            current_symmetry += 1 # Go to next symmetry
        if current_symmetry  == number_symmetries+1:
            break
    return Natural_Occupations
    
    
def Natural_Occupations_CI(load_file, number_symmetries):
    Natural_Occupations = {}
    get_occupations_check = 0 # Check when HF energies are found
    occupations_string = '' # Accumalte HF energies
    occupations_array = np.zeros(1000) # Array to temp store energies
    newline_check = 0 # Check newlines for stop
    current_symmetry = 1
    # CI wavefunctions can have natural occ, for more than one state
    #  need to find the reference state
    ref_state_check = 0 
    for i in load_file:
        if "= the reference state" in i:
            ref_state_check = 1
        if ref_state_check == 1:
            if i[0:11] == " Symmetry "+str(current_symmetry):
                get_occupations_check = 1
                if "No occupied orbitals" in i:
                    get_occupations_check = 0
                    current_symmetry += 1
            elif get_occupations_check == 1:
                if i == "\n":
                    newline_check += 1
                occupations_string = occupations_string + i
            if newline_check == 2:
                get_occupations_check = 0 # Reset
                newline_check = 0 # Reset
                
                data = occupations_string.split(" ")
                counter = 0
                for j in data:
                    if j != '' and j != "\n":
                        occupations_array[counter] = float(j)
                        counter += 1
                Natural_Occupations[current_symmetry] = np.copy(occupations_array[0:counter])
                occupations_string = '' # Reset
                current_symmetry += 1 # Go to next symmetry
            if current_symmetry  == number_symmetries+1:
                break
    return Natural_Occupations
	
	
def metal_d_orbitals(load_file):
    metals = ["Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn",
              "Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
              "Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg",
              "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"]
    orbital_check = False
    metal_d_orbitals = {}
    symmetry_counter = 0
    for line in load_file:
        if "Molecular orbitals for symmetry species" in line:
            symmetry_counter += 1
            metal_d_orbitals[symmetry_counter] = []
        if "    Orbital    " in line and not "":
            orbital_check = True
            # Orbital_number, d1, d2, d3, d4, d5, d_total, total
            orbitals = np.zeros((len(line.split())-1,6))
            for i in range(1, len(line.split())):
                orbitals[i-1,0] = float(line.split()[i])
        elif orbital_check and line == "\n":
            orbital_check = False
            if metal_d_orbitals[symmetry_counter] == []:
                metal_d_orbitals[symmetry_counter] = orbitals
            else:
                metal_d_orbitals[symmetry_counter] = np.vstack((metal_d_orbitals[symmetry_counter], orbitals))
        elif orbital_check:
            line_list = line.split()
            for i in range(3, len(line_list)):
                if line_list[1] in metals and line_list[2] == ":3d2-":
                    orbitals[i-3,1] += abs(float(line_list[i]))
                elif line_list[1] in metals and line_list[2] == ":3d1-":
                    orbitals[i-3,2] += abs(float(line_list[i]))
                elif line_list[1] in metals and line_list[2] == ":3d0":
                    orbitals[i-3,3] += abs(float(line_list[i]))
                elif line_list[1] in metals and line_list[2] == ":3d1+":
                    orbitals[i-3,4] += abs(float(line_list[i]))
                elif line_list[1] in metals and line_list[2] == ":3d2+":
                    orbitals[i-3,5] += abs(float(line_list[i]))
    return metal_d_orbitals
