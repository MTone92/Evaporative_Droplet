# Evaporative_Droplet
This project is designed for calculating the stream function and draw streamlines in droplets induced by evaporation. **Uniform evaporation** and **Nonuniform evaporation** have been coded to calculate the stream function where the uniform mass flux and diffusive mass flux on the droplet surface are considered as boundary conditions, respectively.

## Evaporation_Uniform
The main program for this folder is ```UniformEvporation.m``` which calls other ```*.m``` files as functions and finish the calculation, store the results as a **.mat** file.

Note that this program used parrallel computing in function ```PSI.m``` in which the line 4 ```nodes = 30;``` represents the number of CPU needed is 30. Make sure to change the value of nodes to the proper CPU value on your server.

The ```StreamFunctionPlot.m``` and ```StreamFunctionPlot_Final.m``` are used to plot the streamlines, which is the contour plot of the stream function. They should be executed in order.

To run the code, place the main program and all other ```*.m``` files in one folder.

## Evaporation_Nonuniform
The main program for this folder is ```NonuniformEvaporation.jl``` and all the functions are defined in file ```Module MyFunctions.jl```. Parallel computing is also used in the main program and will allocate all the CPU assigned to this task automatically.

The ```StreamFunctionPlot.m``` and ```StreamFunctionPlot_Final.m``` are used to plot the streamlines, which is the contour plot of the stream function. They should be executed in order.

To run the code, place ```NonuniformEvaporation.jl``` and ```Module MyFunctions.jl``` files in one folder.
