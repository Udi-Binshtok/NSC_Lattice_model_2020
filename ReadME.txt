This .txt file Can be oppend using Matlab software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ReadME file for computer simulation code used in Dray et al. 2020            %%
%% This code was written by Udi Binshtok (udi.binshtok@gmail.com)               %%
%% With acknowledge to Micha Hersch on his contribution to the Morphology parts %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Be sure that you have all of the following in the same folder (in a path that is readable by Matlab):
1) Matlab code file called NSC_lattice_model_simulation.m
2) Matlab code file called CellsPositions.m
3) Matlab code file called NeuronsPositios.m
4) Matlab file called lattice.mat
5) Folder called functions (consist of more Matlab code files and folders)

****************************************************************************************************************************************************************
NOTE: the default number of time steps is 500 (which might take a very long time for the simulation to run, might be 2 days with the default set of parameters).
      Consider changeing it to a lower value (even as low as 5, to start with) in the defaultparams.m file located in the functions folder.  
****************************************************************************************************************************************************************
To run a simulation (with the default parameters) follow these instructions:
1) Open Matlab software.
2) Set Matlab Current Folder to the folder containing the code files and the function folder.
3) In the Command Window type: 
   Load('lattice.mat')
   press enter (to load the file).
4) Choose a Folder for output results and save it as a varible called OutputFolder.
   for example, to save the results into C:\Saved Data Folder\Simulations
   type in the Command Window: 
   OutputFolder = 'C:\Saved Data Folder\Simulations';
   press enter
5) Choose a name for the output result Matlab file and save it as a varaible called OutputName.
   for example, type in the Wommand Window: 
   OutputName = 'simulation_number1';
   press enter
6) To start the simulation type in the Command Window:
   [ geo, Cells, time, Movie, params, Simulation_timing ] = statistic_model_simulation( g, OutputFolder, OutputName, [], [], []);
   press enter

The simulation will start.
Every 10 time steps the current time step will be displayed in the Command Window.
NOTE: If right after the simulation as started there is a problem with a missmatch dimensions error (Matlab error will appear in red in the Command Window), it could be due to a small error with the lattice initial conditions - just try to run the statistic_model_simulation again.

Notes:
1) To change the simulation parameters go to a Matlab code file callled defaultparams.m located in the functions folder.
   In the defaultparams code you can see all of the parameters and change them as you wish.
2) It is possible to provide an outside parameters set and outside initial geometrical and cells structure.
   use [ geo,Cells, time, Movie, params, Simulation_timing ] = statistic_model_simulation( g, OutputFolder, OutputName, params, Initial_geo, Initial_Cells);
   where params is the parameters set (should be similar to the set you get from the function defaultparams.m). 
   Initial_geo and Initial_Cells are the structures geo and Cells (See function NSC_lattice_model_simulation.m for more details).

The positions of the different cell types (including the position of mother cells) during the simulation can be extracted using CellsPositions.m function
The positions of the neurons during the simulation can be extracted using NeuronsPositios.m function

Good luck