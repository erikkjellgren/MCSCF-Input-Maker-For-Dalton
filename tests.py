from Input_Maker import make_input
import numpy as np


def test_MP2_RAS():
    A = make_input.Input_Maker("data/testfiles/mp2_RAS.out")
    A.file_name = "data/from_mp2_RAS.inp"
    A.MCSCF_method = "RAS"
    A.pick_RAS_by_active_threshold(0.01)
    A.write_input_file()
    
    with open("data/testfiles/mp2_RAS_check.inp") as f:
        check_file = list(f)
    with open("data/from_mp2_RAS.inp") as f:
        input_file = list(f)
        
    for i in range(0, len(check_file)):
        assert check_file[i] == input_file[i]
    
    
def test_CI_RAS():
    A = make_input.Input_Maker("data/testfiles/CI_RAS.out")
    A.file_name = "data/from_CI_RAS.inp"
    A.MCSCF_method = "RAS"
    A.pick_RAS_by_active_threshold(0.01)
    A.write_input_file()
    
    with open("data/testfiles/CI_RAS_check.inp") as f:
        check_file = list(f)
    with open("data/from_CI_RAS.inp") as f:
        input_file = list(f)
        
    for i in range(0, len(check_file)):
        assert check_file[i] == input_file[i]
        
        
def test_MP2_RAS2():
    C = make_input.Input_Maker("data/testfiles/mp2_RAS2.out")
    C.MCSCF_method = "RAS"
    C.wavefunction_type = "CI"
    C.file_name = "data/from_CI_RAS2.inp"
    C.inactive = np.array([4, 0])
    C.RAS1 = np.array([8, 2])
    C.RAS2 = np.array([1, 1])
    C.RAS3 = np.array([53, 27])
    C.active_electrons_in_RAS2 = 2
    C.write_input_file()
    
    with open("data/testfiles/CI_RAS2_check.inp") as f:
        check_file = list(f)
    with open("data/from_CI_RAS2.inp") as f:
        input_file = list(f)
        
    for i in range(0, len(check_file)):
        assert check_file[i] == input_file[i]
 
        
def test_lrMCSCF():
    C = make_input.Input_Maker("data/testfiles/lrMP2.out")
    C.MCSCF_method = "CAS"
    C.wavefunction_type = "lrmcscf"
    C.file_name = "data/from_lrMCSCF.inp"
    C.pick_CAS_by_active_threshold(0.01)
    C.write_input_file()
    
    with open("data/testfiles/lrMCSCF.inp") as f:
        check_file = list(f)
    with open("data/from_lrMCSCF.inp") as f:
        input_file = list(f)
        
    for i in range(0, len(check_file)):
        assert check_file[i] == input_file[i]
        
        
def test_occupied_threshold_electron_retrieval():
    mc = make_input.Input_Maker("data/testfiles/Ethene_TZVP.out")
    mc.pick_CAS_occupied_threshold_electron_retrieval(1.9987, retrieval_electron=0.9)
    mc.MCSCF_method="cas"
    mc.write_input_file(check_values_only=True)
    
    CAS_check = np.array([1, 1, 1, 0, 1, 1, 1, 0])
    inactive_check = np.array([2, 0, 1, 0, 2, 0, 0, 0])
    for i in range(0, len(mc.CAS)):
        assert CAS_check[i] == mc.CAS[i]
        assert inactive_check[i] == mc.inactive[i]