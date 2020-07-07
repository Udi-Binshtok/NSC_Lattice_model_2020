function [geo1,Cells1,StopSimulation] = NeuronMigration(geo,Cells,step,params)
%NEURONMIGRATION will migrate Neurons out of the lattice
%   the Neurons at time step 'step' will delaminate out of the lattice at time step 'step' 
   
    StopSimulation = 0;
    
    geo1 = geo;   
    Cells1 = Cells;

    Neurons = find(Cells1(step).states == 5);     % which cells are Neurons at time step 'step'
    
    if ~isempty(Neurons)
        for n = 1:length(Neurons)
            
            g1 = geo1(step).g;
                        
            [ NeuronLocation ] = Location_of_one_cell( g1, Neurons(n) );
            NeuronLocation_x = NeuronLocation(2,2,1);
            NeuronLocation_y = NeuronLocation(2,2,2);
            
            [g1,StopSimulation] = kill_cells(g1,Neurons(n));
            
            if StopSimulation
                return
            end
            Cells1(step).states(Neurons(n)) = -1;   % state of all neurons from time step 'step-1' is '-1' at time step 'step'
            Cells1(step).NeuronMigration(n,1:3) = [ Neurons(n) NeuronLocation_x NeuronLocation_y ]; % record Neuron migration event [ cell_num x_position y_position ]

            geo1(step).g = g1;
            
            
            % find if any cell reduced its bonds to 3, and if so add a bond to this cell 
            for i = 1:(length(g1.cells) - 1)
                g1 = geo1(step).g;
                cell_num = i;
                num_of_bonds = length(g1.cells{cell_num+1});
                if g1.dead(cell_num) == 0 && num_of_bonds < 3 && Cells1(step).states(cell_num) ~= 5
                    disp('weird - the live cell cell_num, which is not a Neuron, has fewer than 3 bonds')
                    StopSimulation = 1;
                    return
%                     keyboard
                end
                if g1.dead(cell_num) == 0 && num_of_bonds == 3
                    [geo1,StopSimulation] = AddBond(geo1,Cells1,step,cell_num);
                    if StopSimulation
                        return
                    end
                end
            end
            
            mi = params.mi.NeuronMigation;
            factor = params.factor;
            geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi,factor);
        end
    end
end