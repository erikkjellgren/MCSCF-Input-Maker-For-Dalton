from Input_Maker import file_loader as fload
from Input_Maker import utility as util
from Input_Maker import analyser as anal
from Input_Maker import heuristics as heu

import numpy as np

class Input_Maker():
    # Bad naming: __type_wavefunction = wavefunction form output file
    #        and: wavefunction_type = wavefunction going into input file
    def __init__(self, output_file):
        """Get numbers from output file"""
        with open(output_file, "r", encoding="utf-8") as f:
            self.__load_file = list(f)
            
        self.__type_wavefunction = fload.wavefunction_type_output(self.__load_file)
        self.__orbitals_in_symmetries = fload.orbital_symmetries(self.__load_file)
        self.__number_symmetreis = len(self.__orbitals_in_symmetries)
        self.__total_nuclei_charge = fload.total_nuclei_charge(self.__load_file)
        
        if self.__type_wavefunction == "ci": # Not MP2 wavefunction
            self.natural_occupations = fload.Natural_Occupations_CI(self.__load_file, self.__number_symmetreis)
            self.__total_electrons = fload.electrons(self.__load_file)
        
        elif self.__type_wavefunction == "mp2":
            self.number_closed_shell = fload.closed_shell_number(self.__load_file)
            self.Hartree_Fock_orbital_energies = fload.HF_orb_energies(self.__load_file, self.__number_symmetreis)
            self.__total_electrons = fload.electronsMP2(self.__load_file)
            self.natural_occupations = fload.Natural_Occupations_MP2(self.__load_file, self.__number_symmetreis)
        
        else:
            self.natural_occupations = fload.Natural_Occupations(self.__load_file, self.__number_symmetreis)
            self.__total_electrons = fload.electrons(self.__load_file)
        
            
        self.natural_occupation_sum = util.Natural_Occupation_Summation(self.natural_occupations)
        self.natural_occupations, self.natural_occupations_index = util.Sort_Natural_Occupations(self.natural_occupations)
        self.relative_natural_occupations_occupied, self.relative_natural_occupations_virtuel = util.Relative_Natural_Occupations(self.natural_occupations)
        
        """Set variables for input files"""
        self.inactive = np.zeros(self.__number_symmetreis, dtype=int)
        self.CAS = np.zeros(self.__number_symmetreis, dtype=int)
        self.RAS1 = np.zeros(self.__number_symmetreis, dtype=int)
        self.RAS2 = np.zeros(self.__number_symmetreis, dtype=int)
        self.RAS3 = np.zeros(self.__number_symmetreis, dtype=int)
        self.RAS1_holes = np.array([0, 2], dtype=int)
        self.RAS3_electrons = np.array([0, 2], dtype=int)
        self.file_name = "input.dal"
        self.spin_multiplicity = 1
        self.MCSCF_method = "undefined"
        self.wavefunction_type = "MCSCF"
        self.symmetry = 1
        self.state = 1
        self.active_electrons = 0
        self.active_electrons_in_RAS2 = 0
        self.max_micro = 24
        self.max_macro = 24
        self.symmetry_threshold = 10**-4
        # MCSCFsrDFT stuff
        self.srxfunctional = "SRXPBEGWS"
        self.srcfunctional = "SRCPBEGWS"
        self.range_separation_parameter = 0.4
        # Response stuff
        self.response = "undifened"
        
        """Some internal variables"""
        self.reorder_neglect_threshold = 0.001
        self.get_nat_occ_neglect_threshold = 0.001
        
        
    def pick_RAS_by_active_threshold(self, threshold):
        self.RAS1, self.RAS3, self.inactive = heu.Pick_RAS_active_threshold(threshold, self.natural_occupations)
        
        
    def pick_CAS_by_active_threshold(self, threshold):
        self.CAS, self.inactive = heu.Pick_CAS_active_threshold(threshold, self.natural_occupations)
        
    def pick_CAS_by_number_occupied(self, number_occupied, allow_more_virtuel=False):
        """
        Function to pick CAS space based on number of wanted occupied active orbitals.
        One virtual orbtial will be added per active orbital.
        
        If 2*number_occupied > number of electrons, more virtual orbitals will still
          be added if allow_more_virtuel=True, allowing for more virtual orbitals than
          occupied orbitals.
          
        Orbitals with occupation number == 2 will never be added.
        """
        self.CAS, self.inactive = heu.Pick_CAS_number_occupied(number_occupied, self.natural_occupations, allow_more_virt=allow_more_virtuel)
        
    def pick_RASCISD_by_number_occupied(self, number_occupied, excitation_from_to=[0,0]):
        """
        Function to pick RASCISD based on number of wanted active occupied orbitals.
        
        Orbitals with occupation number == 2 will never be added.
        """
        self.RAS1, self.RAS2, self.RAS3, self.inactive, self.active_electrons_in_RAS2 = heu.Pick_RASCI_number_occupied(number_occupied, self.natural_occupations, excitation=excitation_from_to)
        
    def scan_threshold_all(self):
        anal.threshold_scan_all(self.natural_occupations)
        
        
    def scan_threshold_per_sym(self):
        anal.threshold_scan_symmetries(self.natural_occupations)
    

    def get_natural_occupancies(self, threshold=2.0):
        """
        Prints the natural occupancies in a formatted manner.
        
        input : threshold, to specif which occupation numbers that should be underlined.
                   will specify with respect to 2. Same number of occupied and 
                   unoccupied will be underlined in each symmetry.
                   default=2.0. I.e. none will be underlined, if 0.0 all will be underlined.
                   
        self.get_nat_occ_neglect_threshold, can be set to decrease or increase the 
        threshold for which occupation numbers will be printed. default=0.001
        """
        anal.print_natural_occ(self.natural_occupations, threshold, self.get_nat_occ_neglect_threshold)
        
        
    def get_relative_natural_occupations(self, show_virtuel_occupations=True):
        anal.print_relative_natural_occ(self.relative_natural_occupations_occupied, self.relative_natural_occupations_virtuel, self.natural_occupations, self.get_nat_occ_neglect_threshold, show_virtuel=show_virtuel_occupations)
        
    def __write_active_space(self):
        self.__input_file.write("*CONFIGURATION INPUT"+"\n")
        self.__input_file.write(".INACTIVE"+"\n")
        self.__input_file.write(" ")
        for i in self.inactive:
            self.__input_file.write(str(i)+" ")
        self.__input_file.write("\n")
        
        if self.MCSCF_method == "cas":
            self.__input_file.write(".CAS SPACE"+"\n")
            self.__input_file.write(" ")
            for i in self.CAS:
                self.__input_file.write(str(i)+" ")
            self.__input_file.write("\n")
        elif self.MCSCF_method == "ras":
            self.__input_file.write(".RAS1 SPACE"+"\n")
            self.__input_file.write(" ")
            for i in self.RAS1:
                self.__input_file.write(str(i)+" ")
            self.__input_file.write("\n")
            
            self.__input_file.write(".RAS2 SPACE"+"\n")
            self.__input_file.write(" ")
            for i in self.RAS2:
                self.__input_file.write(str(i)+" ")
            self.__input_file.write("\n")
            
            self.__input_file.write(".RAS3 SPACE"+"\n")
            self.__input_file.write(" ")
            for i in self.RAS3:
                self.__input_file.write(str(i)+" ")
            self.__input_file.write("\n")
            
            self.__input_file.write(".RAS1 HOLES"+"\n")
            self.__input_file.write(" "+str(self.RAS1_holes[0])+" "+str(self.RAS1_holes[1])+"\n")
            self.__input_file.write(".RAS3 ELECTRONS"+"\n")
            self.__input_file.write(" "+str(self.RAS3_electrons[0])+" "+str(self.RAS3_electrons[1])+"\n")
        
        self.__input_file.write(".ELECTRONS"+"\n")
        self.__input_file.write(" "+str(self.active_electrons)+"\n")
        self.__input_file.write(".SYMMETRY"+"\n")
        self.__input_file.write(" "+str(self.symmetry)+"\n")
        
    def __write_reorder(self):
        symms = [[] for i in range(len(self.natural_occupations))]
        reorder = False
        counter = 0
        for key in self.natural_occupations:
            for i in range(0, len(self.natural_occupations_index[key])):
                if str(i) != str(self.natural_occupations_index[key][i]):
                    if self.natural_occupations[key][i] > self.reorder_neglect_threshold or self.natural_occupations[key][self.natural_occupations_index[key][i]] > self.reorder_neglect_threshold:
                        if self.natural_occupations[key][i] < 2 - self.reorder_neglect_threshold or self.natural_occupations[key][self.natural_occupations_index[key][i]] < 2 - self.reorder_neglect_threshold:
                            symms[counter].append(str(i+1)+" "+str(list(self.natural_occupations_index[key]).index(i)+1))
            counter += 1
        for i in range(len(symms)):
            if len(symms[i]) != 0:
                reorder = True
        
        if reorder == True:
            self.__input_file.write(".REORDER\n")
            for i in range(len(symms)):
                self.__input_file.write(str(len(symms[i]))+" ")
            self.__input_file.write("\n")
            for i in range(len(symms)):
                for j in range(len(symms[i])):
                    self.__input_file.write(str(symms[i][j])+" ")
            self.__input_file.write("\n")
    
    def write_input_file(self):
        """Check coherence in settings"""
        self.MCSCF_method = str(self.MCSCF_method).lower()
        self.wavefunction_type = str(self.wavefunction_type).lower()
        self.response = str(self.response).lower()
        if self.MCSCF_method == "mcscfsrdft":
            self.MCSCF_method = "lrmcscf"
        
        if self.MCSCF_method == "cas":
            self.active_electrons = self.__total_electrons - np.sum(2*self.inactive, dtype=int) # number of active electrons must be this. Cannot infer number electrons just from CAS
        elif self.MCSCF_method == "ras":
            self.active_electrons = np.sum(2*self.RAS1, dtype=int) + self.active_electrons_in_RAS2 # Cannot infer from output file how many active electrons should be accounted for in RAS2
        else:
            assert False, "MCSCF_method is invalid, choose RAS or CAS"
        
        assert self.active_electrons ==  self.__total_electrons - np.sum(2*self.inactive, dtype=int), "Number of active electrons does not match active/inactive space"
        assert self.active_electrons%2 == 0, "Number of active electrons is uneven"
        assert self.__total_electrons%2 == 0, "Number of total electrons is uneven"
        
        """Print something useful"""
        print("Total molecular charge:", self.__total_nuclei_charge - self.__total_electrons)
        
        """Write input file"""
        self.__input_file = open(self.file_name, "w+")
        self.__input_file.write("**DALTON INPUT"+"\n")
        if self.response == "excitation" or self.response == "excitation_tda":
            self.__input_file.write(".RUN RESPONSE\n")
        else:
            self.__input_file.write(".RUN WAVE FUNCTION"+"\n")
        self.__input_file.write("*MOLBAS\n")
        self.__input_file.write(".SYMTHR\n")
        self.__input_file.write(" "+str(self.symmetry_threshold)+"\n")
        
        """Write CI"""
        if self.wavefunction_type == "ci":
            self.__input_file.write("**WAVEFUNCTION"+"\n")
            self.__input_file.write(".CI\n")
            self.__input_file.write("*CI VECTOR"+"\n")
            self.__input_file.write(".PLUS COMBINATIONS"+"\n")
            self.__input_file.write("*CI INPUT"+"\n")
            self.__input_file.write(".CINO"+"\n")
            self.__input_file.write(".MAX ITERATIONS\n")
            self.__input_file.write(" 100\n")
            self.__input_file.write(".STATE"+"\n")
            self.__input_file.write(" "+str(self.state)+"\n")
            
            self.__write_active_space()
            
            self.__input_file.write("*ORBITAL INPUT"+"\n")
            self.__input_file.write(".MOSTART"+"\n")
            self.__input_file.write(" NEWORB"+"\n")
            
            self.__write_reorder()
            
            self.__input_file.write("*OPTIMIZATION"+"\n")
            self.__input_file.write(".DETERMI"+"\n")
            self.__input_file.write(".MAX MICRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_micro)+"\n")
            self.__input_file.write(".MAX MACRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_macro)+"\n")
            
        
        """Write MCSCF"""
        if self.wavefunction_type == "mcscf":
            self.__input_file.write("**WAVEFUNCTION"+"\n")
            self.__input_file.write(".MCSCF\n")
            
            self.__write_active_space()
            
            self.__input_file.write("*ORBITAL INPUT"+"\n")
            self.__input_file.write(".MOSTART"+"\n")
            self.__input_file.write(" NEWORB"+"\n")
            
            self.__write_reorder()
            
            self.__input_file.write("*OPTIMIZATION"+"\n")
            self.__input_file.write(".STATE"+"\n")
            self.__input_file.write(" "+str(self.state)+"\n")
            self.__input_file.write(".MAX MICRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_micro)+"\n")
            self.__input_file.write(".MAX MACRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_macro)+"\n")
        
        
        """Write lr-MCSCF"""
        if self.wavefunction_type == "lrmcscf":
            self.__input_file.write("**WAVEFUNCTION"+"\n")
            self.__input_file.write(".MCSRDFT\n")
            self.__input_file.write(".SRFUN\n")
            self.__input_file.write(" "+str(self.srxfunctional)+" "+str(self.srcfunctional)+"\n")
            
            self.__write_active_space()
            
            self.__input_file.write("*ORBITAL INPUT"+"\n")
            self.__input_file.write(".MOSTART"+"\n")
            self.__input_file.write(" NEWORB"+"\n")
            
            self.__write_reorder()
            
            self.__input_file.write("*OPTIMIZATION"+"\n")
            self.__input_file.write(".STATE"+"\n")
            self.__input_file.write(" "+str(self.state)+"\n")
            self.__input_file.write(".MAX MICRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_micro)+"\n")
            self.__input_file.write(".MAX MACRO ITERATIONS\n")
            self.__input_file.write(" "+str(self.max_macro)+"\n")
            self.__input_file.write("**INTEGRALS\n")
            self.__input_file.write("*TWOINT\n")
            self.__input_file.write(".DOSRIN\n")
            self.__input_file.write(".ERF\n")
            self.__input_file.write(" "+str(self.range_separation_parameter)+"\n")
            
        """Write singlet excitation response"""
        if self.response == "excitation" or self.response == "excitation_tda":
            self.__input_file.write("**RESPONSE\n")
            if self.response == "excitation_tda":
                self.__input_file.write(".TDA\n")
            self.__input_file.write("*LINEAR\n")
            self.__input_file.write(".SINGLE RESIDUE\n")
            self.__input_file.write(".ROOTS\n")
            for i in range(self.__number_symmetreis):
                self.__input_file.write(" 3")
            self.__input_file.write("\n")
            
            self.__input_file.write(".NSTART\n")
            for i in range(self.__number_symmetreis):
                self.__input_file.write(" 10")
            self.__input_file.write("\n")
            
            self.__input_file.write(".NSIMUL\n")
            for i in range(self.__number_symmetreis):
                self.__input_file.write(" 6")
            self.__input_file.write("\n")
            self.__input_file.write(".PRINT\n")
            self.__input_file.write(" 4\n")
            
            
        self.__input_file.write("**END OF DALTON INPUT"+"\n")
        self.__input_file.close()
        
       
