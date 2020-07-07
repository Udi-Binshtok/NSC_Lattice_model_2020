function [geo,Cells,time,StopSimulation,Simulation_timing] = simulation(Initial_geo,Initial_Cells,params)

%   Each cell is at one of these states:
%   NSC (0), aNSC (1), MotherCell (2), DaughterCell (3), aNP (4), div_aNP (4.1), Neuron (5), deleted cell/migrated Neuron (-1)
%   state 4.1 is a div_aNP cell which is an aNP that has divided once and will differenitate on the next step (it will be presented as aNP)
   
    T = params.Time_max;               
    time = (1:1:T);
    
    %%% These are the initial time steps
    Cells = Initial_Cells;
    geo = Initial_geo;
    
    %%% make sure that all vertices have only 3 sharing cells 
    pert = 0.20;
    [geo,StopSimulation] = ReduceVertexSharingCellNumber(geo,Cells,1,pert);
    if StopSimulation
        return
    end
    mi = params.mi.simulation;
    factor = params.factor;
    geo = MorphologyUpdate(geo,Cells,1,params,mi,factor);
        
    %%% simulation loop      
    for step = 2:1:T
        
        %%% display simulation progress  
        if mod(step,10) == 0 % display the current time step every 10 time steps 
            disp(['time step ' num2str(step)])
        end
        
        %%% setting default geo.g and Cells.states at time step 'step', and update it in the next set of functions
        geo(step).g = geo(step-1).g; 
        Cells(step).states = Cells(step-1).states;
        
        %%% NSC to aNSC (regulated by neighbors aNP) or directly to Neuron
        nsc = tic;
        [geo,Cells] = NSC2aNSCor2Neuron(geo,Cells,step,params);
Simulation_timing(step).NSC = toc(nsc);
        
        %%% aNSC to Mother Cell
        ansc = tic;
        [geo,Cells] = aNSC2MotherCell(geo,Cells,step,params);
Simulation_timing(step).aNSC = toc(ansc);

        %%% Mother cells to Daughter cells
        mc = tic;
        [geo,Cells,StopSimulation] = MotherCells2DaughterCells(geo,Cells,step,params);
        if StopSimulation
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            return
        end
        for i = 1:(length(geo(step).g.cells) - 1) 
            num_of_bonds(i) = length(geo(step).g.cells{i+1});
        end
        if ~isempty(find(num_of_bonds == 2))
            disp('problem - some cells have only 2 bonds')
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            StopSimulation = 1;
            return
        end
Simulation_timing(step).MC = toc(mc);        
                    
        %%% Daughter cells to NSC or aNP 
        DC = tic;
        [geo,Cells] = DaughterCellsFate(geo,Cells,step,params);
Simulation_timing(step).DC = toc(DC);
        
        %%% aNP to Neurons or aNP division 
        anp = tic;
        [geo,Cells,StopSimulation] = aNPsFate(geo,Cells,step,params);
        if StopSimulation
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            return
        end
        for i = 1:(length(geo(step).g.cells) - 1) 
            num_of_bonds(i) = length(geo(step).g.cells{i+1});
        end
        if ~isempty(find(num_of_bonds == 2))
            disp('problem - some cells have only 2 bonds')
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            StopSimulation = 1;
            return
        end
Simulation_timing(step).aNP = toc(anp);
                    
        %%% Neuron migration
        neur = tic;
        [geo,Cells,StopSimulation] = NeuronMigration(geo,Cells,step,params);
        if StopSimulation
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            return
        end
        for i = 1:(length(geo(step).g.cells) - 1) 
            num_of_bonds(i) = length(geo(step).g.cells{i+1});
        end
        if ~isempty(find(num_of_bonds == 2))
            disp('problem - some cells have only 2 bonds')
            if ~exist('Simulation_timing')
                Simulation_timing = NaN;
            end
            StopSimulation = 1;
            return
        end
Simulation_timing(step).Neuron = toc(neur);
        
        % make sure that all vertices has only 3 sharing cells
        pert = 0.20;
        [geo,StopSimulation] = ReduceVertexSharingCellNumber(geo,Cells,step,pert);
        if StopSimulation
            return
        end
        morph = tic;
        mi = params.mi.simulation;
        factor = params.factor;
        geo = MorphologyUpdate(geo,Cells,step,params,mi,factor);
Simulation_timing(step).morphology = toc(morph);
        
    end
        
end