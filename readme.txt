File description:
1. Model.m: matlab function including stochastic simulations of Msn2 and Msn4 dynamics based on dynamics of PKA and X, identification of Msn2/4 pulses and calculation of pulse coincidence rate.
2. Simulation.m: matlab code for simulating cells under different parameters and calculating pulse coincidence rate based on Model.m. (Expected run time for Simulation.m on a "normal" desktop computer is 4 h, note that cell number for simulation is 3. For more cells, parallel running is recommended)
3. expected output: example output of Simulation.m. In each subfolder, csv file includes the average results of these 3 cells and pdf files are Msn2/Msn4/PKA/X trajectories of these 3 cells.

System requirements:
1. Matlab of most versions (e.g., R2020a in Windows10 x64 system)

Installation guide: (Typical install time on a "normal" desktop computer: a few seconds)
1. Add path of Model.m to Matlab (as shown in Simulation.m)