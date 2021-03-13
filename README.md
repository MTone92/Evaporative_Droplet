# Evaporative_Droplet
This project is designed for calculating the stream function and draw streamlines in droplets induced by evaporation. **Uniform evaporation** and **Nonuniform evaporation** have been coded to calculate the stream function where the uniform mass flux and diffusive mass flux on the droplet surface are considered as boundary conditions, respectively.
## Evaporation_Uniform
The main program for this folder is *UniformEvporation.m* which calls other *.m* files as functions and finish the calculation, store the results as a *.mat* file. Note that this program used parrallel computing in function *PSI.m* in which the line 4 ```nodes = 30;```
