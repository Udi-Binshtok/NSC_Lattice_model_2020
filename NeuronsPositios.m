function [Neurons_positions_StructureArray] = NeuronsPositios(FilesFolder,generic_name,input_simulations,min_time_step,max_time_step)
%NEURONSPOSITIONS will find the positions of neurons in a 3D volume consist of layers of lattices from time step min_time_step to time step max_time_step 
%   input:
%   FilesFolder is the folder where the saved simulations data (from the statistic_model_simulation function) is kept as a .mat file (example: FilesFolder = 'C:\Saved Data Folder\Simulations';) 
%       (see also OutputFolder in statistic_model_simulation)
%   generic_name is a name that is shared to all of the simulations files, in case that you generated more than one simulation (example: generic_name = 'simulation_number') 
%       (see also OutputName in statistic_model_simulation)
%   input_simulations is the number of the simulations that you want this code to go through (example: input_simulations = 1:5; % will read simulation_number1, simulation_number2, ..., simulation_number5)  
%   min_time_step is the lowest time step from which the code will find the neurons positions (example: min_time_step = 135;) 
%   max_time_step is the highest time step from which the code will find the neurons positions (example: max_time_step = 500;) 
%   
%   output:
%   Neurons_positions_StructureArray is a structure array where each line is a simulation (from the input simulations) that consist of array of the neurons positions
%       The neurons positions array: columns 1, 2 and 3 are the x,y and z positions (in microns)  

    for s = input_simulations
        FilePath = [FilesFolder '\' generic_name num2str(s) '.mat'];
        input = load(FilePath);
        Cells = input.Cells;
        geo = input.geo;
        
        Z_position_inDays = 1;
        Z_position_inMicrons = Z_position_inDays*0.5; % neurons decent 30 microns in 2 months 
        Neurons_positions = [];
        
        for time_step = min_time_step:max_time_step
            Neurons = Cells(time_step).NeuronMigration;
            if ~isempty(Neurons)
                g = geo(time_step).g;
                % estimate cell diameter at time_step
                cell_area = zeros(1,length(g.dead(g.dead == 0)));
                not_dead = find(g.dead == 0);
                for c = 1:length(not_dead)
                    cell_area(c) = cellarea(g,not_dead(c));
                end
                average_cell_diameter = 2*sqrt(mean(cell_area./pi));
                
                Neurons_positions_xy = Neurons(:,2:3);
                Neurons_positions_xy_inCellDiameter = Neurons_positions_xy./average_cell_diameter;
                Neurons_positions_xy_inMicrons = 10.*Neurons_positions_xy_inCellDiameter; % 1 cell diameter is about 10 microns.
                Neurons_positions_z_inMicrons = Z_position_inMicrons*ones(size(Neurons_positions_xy_inMicrons,1),1);
                Neurons_positions = [ Neurons_positions ; Neurons_positions_xy_inMicrons Neurons_positions_z_inMicrons ];
            end
            Z_position_inDays = Z_position_inDays + 1;
            Z_position_inMicrons = Z_position_inDays*0.5;
        end
        
        Neurons_positions_StructureArray(s).Neurons_positions_xyz = Neurons_positions;
        
    end
end

