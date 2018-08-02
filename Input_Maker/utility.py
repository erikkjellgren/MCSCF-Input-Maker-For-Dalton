import numpy as np


def Natural_Occupation_Summation(Natural_Occupations):
    nat_occ_sum = np.zeros(len(Natural_Occupations))
    counter = 0
    for key in Natural_Occupations:
        nat_occ_sum[counter] = np.sum(Natural_Occupations[key])
        counter += 1
    return nat_occ_sum
    
    
def Sort_Natural_Occupations(Natural_Occupations_in):
    natural_occupations_index = {}
    natural_occupations = {}
    for key in Natural_Occupations_in:
        natural_occupations[key] = np.sort(Natural_Occupations_in[key])[::-1]
        natural_occupations_index[key] = np.argsort(Natural_Occupations_in[key])[::-1]
        
    for key in Natural_Occupations_in: # Nat occ == 2.0, would get unorded index from argsort
        for i in range(0, len(Natural_Occupations_in[key])):
            if Natural_Occupations_in[key][i] != 2.0:
                break
            natural_occupations_index[key][i] = i
    return natural_occupations, natural_occupations_index
    
    
def Relative_Natural_Occupations(Natural_Occupations):
    """
    Takes in a dictionary of natural occupations, and returns a dictionary of
    relative occupied natural occupations, and a dictionary of relative 
    virtuel natural occupations.
    
    n_rel_occ = (2 - n_occ)/(2 - min(n_occ))
    n_rel_virt = n_virt/max(n_virt)
    """
    relative_occupied = {}
    relative_virtuel = {}
    rel_occ = 0
    rel_virt = 0
    idx_counter = 0
    for key in Natural_Occupations:
        for occupation in Natural_Occupations[key]:
            if occupation > 1.0 and 2 - occupation > rel_occ:
                rel_occ = 2 - occupation
            elif occupation < 1.0 and occupation > rel_virt:
                rel_virt = occupation
        idx_counter += 1
    idx_counter = 0
    for key in Natural_Occupations:
        relative_occupied[key] = []
        relative_virtuel[key] = []
        for occupation in Natural_Occupations[key]:
            if occupation > 1.0:
                relative_occupied[key].append((2-occupation)/rel_occ)
            else:
                relative_virtuel[key].append(occupation/rel_virt)
        idx_counter += 1
    return relative_occupied, relative_virtuel
    
    
    
    
    
    
    
    