function [geo1, Cells1,StopSimulation] = MotherCells2DaughterCells(geo, Cells, step, params)
%MOTHERCELLS2DAUGHTERCELLS will turn Mother cells at time step 'step-1' to Daughter cells at time step 'step'
%   Mtoher cells will divide into 2 Daughter cells 
   
    StopSimulation = 0;
    
    geo1 = geo;
    g = geo1(step).g;
    Cells1 = Cells;

    MotherCells = find(Cells1(step).states == 2); % which cells are MotherCells at time step 'step'
    if ~isempty(MotherCells)
        new_cells_num = NaN(1,length(MotherCells));
        for m = 1:length(MotherCells)
            [ MCLocation ] = Location_of_one_cell( g, MotherCells(m) );
            MCLocation_x = MCLocation(2,2,1);
            MCLocation_y = MCLocation(2,2,2);
            [g, new_cells_num(m), weird] = cell_division( g, MotherCells(m) ); % divide Mother cells
            if weird
                disp('weird - some cells have only 2 bonds')
                StopSimulation = 1;
                return
%                 keyboard
            end
        Cells1(step).states(MotherCells(m)) = 3;   % state of all Mother cells from time step 'step-1' is '3' at time step 'step'
        Cells1(step).states(new_cells_num(m)) = 3; % state of all new cells is '3'
        Cells1(step).MC_division(m,1:4) = [ MotherCells(m), new_cells_num(m), MCLocation_x, MCLocation_y ]; % record Mother cell (MC) division into 2 Daughter cells (DC) [ DC1, DC2, MC_position_x, MC_position_y ] ] 
        g.dead(new_cells_num(m)) = 0;              % new cells are not dead
        g.xboundary(new_cells_num(m),1:2) = 0;     % 0 means periodic boundary condition
        g.yboundary(new_cells_num(m),1:2) = 0;     % 0 means periodic boundary condition
        geo1(step).g = g;
        
        % special case morphology update - increasing daugher cell size 
        mi_sp = params.mi.MotherCells2DaughterCells_sp;
        factor = params.factor;
        geo1 = special_MorphologyUpdate(geo1,Cells1,step,params,mi_sp,[MotherCells(m) new_cells_num(m)],factor);
        
        mi = params.mi.MotherCells2DaughterCells;
        factor = params.factor;
        geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi,factor);
        end
        
    end
    
end

