function [] = CellsPositions(FilesFolder,OutputFolder,generic_name,input_simulations,range_steps)
%CELLSPOSITIONS will find and save the positions of cells (from different types, including mother cell (MC) positions) at each time step within the range_steps 
%   input:
%   FilesFolder is the folder where the saved simulations data (from the statistic_model_simulation function) is kept as a .mat file (example: FilesFolder = 'C:\Saved Data Folder\Simulations';) 
%       (see also OutputFolder in statistic_model_simulation)
%   OutputFolder is the folder you wish to save the positions structure array into example: OutputFolder = 'C:\Saved Data Folder\positions'; 
%   generic_name is a name that is shared to all of the simulations files, in case that you generated more than one simulation (example: generic_name = 'simulation_number') 
%       (see also OutputName in statistic_model_simulation)
%   input_simulations is the number of the simulations that you want this code to go through (example: input_simulations = 1:5; % will read simulation_number1, simulation_number2, ..., simulation_number5)  
%   range_steps is the time steps you wish to get the cells positions from (example: range_steps = 300:500;) 


    AddFunctionsFolder = genpath('functions'); % The functions folder contain all of the functions to be used in the simulation  
    addpath(AddFunctionsFolder)
    
    for s = input_simulations
        
        FilePath = [FilesFolder '\' generic_name num2str(s) '.mat'];
        input = load(FilePath);
        Cells = input.Cells;
        geo = input.geo;
        
        for t = range_steps(1):range_steps(end)

            time_step = t;

            qNSC = find(Cells(time_step).states == 0);
            aNSC = find(Cells(time_step).states == 1);
            MC = Cells(time_step).MC_division;
            if ~isempty(MC)
                MC = Cells(time_step).MC_division(:,1)';
            end
            aNP = find(Cells(time_step).states == 4);
            div_aNP = find(Cells(time_step).states == 4.1);

            g = geo(time_step).g;

            Positions(time_step).qNSC = NaN(size(qNSC,2),2); % Each row is a cell. Column 1 is x coordinate, column 2 is y coordinate 
            Positions(time_step).aNSC = NaN(size(aNSC,2),2);
            Positions(time_step).MC = NaN(size(MC,2),2);
            Positions(time_step).aNP = NaN(size(aNP,2),2);
            Positions(time_step).div_aNP = NaN(size(div_aNP,2),2);

            % filling in qNSC positions:
            for i = 1:size(qNSC,2)
                cell_number = qNSC(i);
                [ cell_periodic_boundary_condition_location ] = Location_of_one_cell( g, cell_number );
                cell_location = cell_periodic_boundary_condition_location(2,2,:);
                Positions(time_step).qNSC(i,1:2) = cell_location(1:2);
            end

            % filling in aNSC positions:
            for i = 1:size(aNSC,2)
                cell_number = aNSC(i);
                [ cell_periodic_boundary_condition_location ] = Location_of_one_cell( g, cell_number );
                cell_location = cell_periodic_boundary_condition_location(2,2,:);
                Positions(time_step).aNSC(i,1:2) = cell_location(1:2);
            end

            % filling in MC positions:
            for i = 1:size(MC,2)
                Positions(time_step).MC(i,1:2) = Cells(time_step).MC_division(i,3:4);
            end

            % filling in aNP positions:
            for i = 1:size(aNP,2)
                cell_number = aNP(i);
                [ cell_periodic_boundary_condition_location ] = Location_of_one_cell( g, cell_number );
                cell_location = cell_periodic_boundary_condition_location(2,2,:);
                Positions(time_step).aNP(i,1:2) = cell_location(1:2);
            end

            % filling in div_aNP positions:
            for i = 1:size(div_aNP,2)
                cell_number = div_aNP(i);
                [ cell_periodic_boundary_condition_location ] = Location_of_one_cell( g, cell_number );
                cell_location = cell_periodic_boundary_condition_location(2,2,:);
                Positions(time_step).div_aNP(i,1:2) = cell_location(1:2);
            end
        end
        
        OutputName = ['Positions_simulation' num2str(s)];    
        if ~isempty(OutputName)
            SaveTo = [OutputFolder '\' OutputName '.mat'];
            save(SaveTo,'Positions')
        end
        
    end
    
end