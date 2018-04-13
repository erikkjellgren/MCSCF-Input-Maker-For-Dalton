# MCSCF Input Maker For Dalton

Script to help making and checking sanity of Dalton input files, for CASSCF or RASSCF methods.

## Example:

As an input the script takes the output file from a previous Dalton calculation. 

```python
C = make_input.Input_Maker("some_output_file.out")
C.MCSCF_method = "CAS"
C.wavefunction_type = "MCSCF"
C.pick_CAS_by_active_threshold(0.01)
C.file_name = "test.dal"
C.write_input_file()
>>>
Total molecular charge: 0
```

## Inspect Natural orbitals:

```python
C.scan_threshold_per_sym()
>>>
Symmetry:  1
Threshold: 0.0220 gives:  1 occupied-active orbitals
Threshold: 0.0194 gives:  2 occupied-active orbitals
Threshold: 0.0183 gives:  3 occupied-active orbitals
Threshold: 0.0171 gives:  4 occupied-active orbitals
Threshold: 0.0161 gives:  5 occupied-active orbitals
Threshold: 0.0135 gives:  6 occupied-active orbitals
Threshold: 0.0126 gives:  7 occupied-active orbitals
Threshold: 0.0095 gives:  8 occupied-active orbitals
Threshold: 0.0047 gives:  9 occupied-active orbitals
Symmetry:  2
Threshold: 0.0415 gives:  1 occupied-active orbitals
Threshold: 0.0317 gives:  2 occupied-active orbitals
```

```python
C.get_natural_occupancies(0.01)
>>>
Symmetry:  1

1.9989   1.9979   1.9979   1.9976   1.9973   
1.9970   1.9964   1.9959   1.9959   1.9951   
1.9947   1.9945   1.9943   1.9938   1.9923   
1.9922   0.0060   0.0053   0.0047   0.0042   
0.0039   0.0039   0.0037   0.0027   0.0024   
0.0022   0.0022   0.0021   0.0021   0.0020   
0.0019   0.0019   0.0017   0.0017   0.0015   

Symmetry:  2

1.9973   1.9949   1.9908   1.9846   1.9737   
1.9728   0.0246   0.0238   0.0124   0.0063   
0.0024   0.0023   0.0021   0.0017   0.0015   
```
