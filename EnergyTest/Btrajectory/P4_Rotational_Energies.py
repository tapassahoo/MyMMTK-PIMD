#  SUMMARY

#  These are values outputted by the rotational energy estimator in one of the test-pimd-energy_analysis.py file's print statements. 

#  This test-pimd-energy_analysis.py file imports a .nc file generated by our main simulation script, and these values were
#  created using 4 beads (P = 4) and 100 000 steps for the first column and 4 beads (P = 4) and 1 000 000 steps for the 
#  second column. The third column consists of the exact energies calculated in the average_energy_test2.py file.

#  The system was 1 translationless H2O rigit rotor, in this set of simulations tau decreased by the amount of beads was kept constant.

#  The simuation was repeated for 4 temperatrues, T = [10, 20, 30, 50] and the data from the simulations is shown above.
 

#  RESULTS 
        
# MoRiBS Results:   Temp    Erot    Error || MMTK 2e5 steps:   Temp    Erot      Error || MMTK 2e6 steps:   Temp    Erot      Error || MMTK 2e6 steps, 150 skipped steps:   Temp    Erot      Error 
                                                                                                                                                                                                                                                              
                    10.0   4.338763                            10.0   4.24513   0.04974                     10.0   4.1465    0.0157                                         10.0   4.5134    0.0159
                    20.0   22.88446                            20.0   23.0997   0.09999                     20.0   22.638    0.0317                                         20.0   23.732    0.0305
                    30.0   39.46542                            30.0   39.9204   0.14435                     30.0   39.894    0.0454                                         30.0   38.611    0.0459
                    50.0   71.99282                            50.0   70.2389   0.24125                     50.0   69.544    0.0752                                         50.0   72.758    0.0709
                                                                                                     