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
    