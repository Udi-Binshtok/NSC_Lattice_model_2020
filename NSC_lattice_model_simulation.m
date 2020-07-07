function [ geo,Cells, time, Movie, params, Simulation_timing ] = NSC_lattice_model_simulation( g, OutputFolder, OutputName, params, Initial_geo, Initial_Cells)
%NSC_LATTICE_MODEL_SIMULATION is the main code to run a computer simulation for Dray et al. 2020 
%   This computer simulation code was written by Udi Binshtok (udi.binshtok@gmail.com)
%   Please first read the ReadMe.txt file for more details.
%   
%   Input: 
%   g is the geometric properties structure of the 2D vertex model lattice of cells (g is stored in the lattice.mat file) 
%   OutputFolder is the fullpath for a folder for saved simulation data (example: OutputFolder = 'C:\Saved Data Folder\Simulations';)  
%   OutputName is the name of the file that the data will be saved under, in the Output folder (example: OutputName = 'simulation_number1';) 
%   
%   Optional input (1. params  2. Initial_geo and Initial_Cells (must use them together): 
%   Note: if the one or more of the optional inputs is not provided repalce it with [] (empty array) 
%   you can run the code as Laure_statistic_model_disordered( g, OutputFolder, OutputName, [], [], [])):
%   params is a provided parameters set fro the user (must contain the values for all the parameters as in defaultparams.m) 
%   Initial_geo is an initial condition for geo, which replaces the Initial_geo structure from get_IC.m function 
%   Initial_Cells is an initial condition for Cells, which replaces the Initial_Cells structure from get_IC.m function
%   
%   Output:
%   geo is the geometric properties structure, such that - 
%       geo(t).g is the geometric properties structure of the lattice of cells at time step 't'
%   Cells is Cells states and events information such that - 
%       Cells(t).states is a vector consist of the state of each cell at time step 't' (including cells that has already delaminated from the lattice) 
%       Cells(t).NSC2Neuron is a vector consist of the qNSC (cell number) that directly differentiate into neurons at time step 't'
%       Cells(t).MC_division is an array consist of divisions events of aNSC at time step 't' where column 1 is the aNSC number that has divided, column 2 is the number of the new cell, columns 3 and 4 are the position x and y of the aNSC that has divided    
%       Cells(t).NeuronMigration is an array consist differentiation events of qNSCs and aNPs at time step 't' where column 1 is the cell number that has differentiated, columns 2 and 3 are the position x and y of the cell that has differentiated 
%   time is a vector consist of the time steps in the simulation
%   Movie is a structure consist of the movie data
%   params is a structure consist of the parameters that were used in the simulation 
%   Simulation_timing is a structure with timing of different parts of the code 
%   
%   During the simulation, at time step t, each cell is characterized by its geometry in geo(t).g and by its state in Cells(t).states 
%   Each cell is at one of these states:
%   NSC (0), aNSC (1), MotherCell (2), DaughterCell (3), aNP (4), div_aNP (4.1), Neuron (5), deleted cell/migrated Neuron (-1)
%   state 4.1 is a div_aNP cell which is an aNP that has divided once and will differenitate on the next step (it will be presented as aNP)
    
    plotMovie = 1;
%     plotMovie = 0; uncomment if you wish not to play the simulated movie 
    
    StopSimulation = 0; % in case of any problem during the simulation StopSimulation would turn to 1 and the simulation will stop with an error message (otherwise it is 0 and the simulation will keep running until it is done)  
    
    StartTiming = tic; % simulation timing 

    AddFunctionsFolder = genpath('functions'); % The functions folder contain all of the functions to be used in the simulation  
    addpath(AddFunctionsFolder)


    %%% Initial conditions
    if ~isempty(Initial_geo) && ~isempty(Initial_Cells)
        %%% parameters
        if isempty(params) 
            [ params ] = defaultparams(Initial_geo.g); % get the default parameters if none provided
        end
    else
        %%% parameters
        if isempty(params)
            [ params ] = defaultparams(g); % get the default parameters if none provided
        end
        mi = params.mi;
        [Initial_geo,Initial_Cells,StopSimulation] = get_IC(g,params,mi.get_IC);
    end
    
    Initial_geo(1).g.paras = params.g_paras; % set the factors for the mechanical energy minimization 
    
    if StopSimulation == 0    
        %%% simulation
        [geo,Cells,time,StopSimulation,Simulation_timing] = simulation(Initial_geo,Initial_Cells,params);
        if StopSimulation
            disp([OutputName ' has stopped'])
        else
            disp([OutputName ' is done']);
        end
    
    else
        Cells = Initial_Cells;
        geo = Initial_geo;
        T = params.Time_max;            
        time_step = params.time_step;    
        time = (1:time_step:T);
        disp([OutputName ' has stopped at the get_IC stage'])
    end
    
SimulationTime = toc(StartTiming);
disp(['The simulation running time : ' num2str(SimulationTime) ' seconds'])
    
    %%% save the simulated data
    if ~isempty(OutputName)
        SaveTo = [OutputFolder '\' OutputName '.mat'];
        save(SaveTo,'geo','Cells','time','params','Simulation_timing')
    end
    
    %%% plot
    if plotMovie
        [Movie,MovieName] = plotting(geo,Cells,time);
        implay([MovieName '.avi'])
    else
        Movie = NaN;
    end
    
end
