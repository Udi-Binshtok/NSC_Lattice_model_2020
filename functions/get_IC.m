function [Initial_geo,Initial_Cells,StopSimulation] = get_IC(g,params,mi)
%GET_IC sets the initial geometry, geo, and initial cells' state, Cells.states
%   Each cell is at one of these states:
%   NSC (0), aNSC (1), MotherCell (2),DaughterCell (3), aNP (4), div_aNP (4.1), Neuron (5), deleted cell/migrated Neuron (-1)
    
    StopSimulation = 0;

    per.aNSC = params.per_aNSC; % observed percentage of aNSC in the Medial pallium (DM)
    per.aNP = params.per_aNP;   % observed percentage of aNP in the Medial pallium (DM)
    
    %%% initialy set all the cells to be NSC (and change it afterwards)
    k = length(g.cells) - 1;            % number of cells in the lattice + dead cells
    Cells(1).states = zeros(1,k);  % set all cells to be at state 0 (NSC)
    
    %%% set all dead cells state to be -1 (deleted cell/migrated Neuron)
    Cells(1).states(g.dead == 1) = -1;
    num_of_cells_in_the_lattice = length(find(~g.dead));
    
    %%% insert aNSCs
    num_of_aNSC = floor((per.aNSC/100)*num_of_cells_in_the_lattice);     
    r = pi*rand(1,num_of_aNSC);                 % get random radi for the position of each aNSC (min r = 0, max r = pi)
    teta = linspace(0,2*pi,num_of_aNSC + 1);    % set the aNSC cells apart from each other in different angle region
    pos = NaN(2,num_of_cells_in_the_lattice);   % get the position of all of the cells in the lattice (first row is r, second row is teta)
    cells_in_lattice = find(~g.dead);
    for c = 1:num_of_cells_in_the_lattice 
        [ full_position ] = Location_of_one_cell( g, cells_in_lattice(c) ); % position of cell c in the lattice (and also in periodic_boundary_condition lattices)
        pos(1,c) = full_position(2,2,3);    % this is r of cell c in the lattic
        pos(2,c) = full_position(2,2,4);    % this is angle of cell c in the lattic
    end
    negative_angle = find(pos(2,:) < 0);
    for t = negative_angle
        pos(2,t) = pi + (pi + pos(2,t));    % map negative angles to positive, that is (0,-pi) -> (2pi,pi)
    end
    aNSC_cell_num = NaN(1,num_of_aNSC);
    for n = 1:num_of_aNSC  % choose the aNSC cells to be closest to r at each of teta regions (one aNSC at each region)                 
        in_teta_region = find(pos(2,:) > teta(n) & pos(2,:) < teta(n+1));
        r_difference = abs(pos(1,:) - r(n).*ones(1,num_of_cells_in_the_lattice));
        new_aNSC = in_teta_region(r_difference(in_teta_region) == min(r_difference(in_teta_region)));
        if length(new_aNSC) > 1
            new_aNSC = new_aNSC(1);
        end
        aNSC_cell_num(n) = cells_in_lattice(new_aNSC);
    end
    Cells(1).states(aNSC_cell_num) = 1;    
    
    %%% insert aNP cells (cannot neighbor with aNSC)  
    % neighbors of aNSC cells
    aNSC = find(Cells(1).states == 1);
    Neighbors = [];
    for n = 1:num_of_aNSC
        neighbors = g.bonds(g.cells{aNSC(n)+1},4)';     % neighbors of one aNSC
        Neighbors = [Neighbors neighbors];              % all of the neighbors of all of the aNSC
    end
    % candidates for aNP
    possible_cells = find(Cells(1).states == 0);
    possible_cells = possible_cells(~ismember(possible_cells,Neighbors));    % cannot neighbor with aNSC
    % set aNP
    num_of_aNP = floor((per.aNP/100)*num_of_cells_in_the_lattice);
    num_of_possible_cells = length(possible_cells);
    random = 1 + (num_of_possible_cells-1)*rand(1,num_of_aNP); 
    Cells(1).states(possible_cells(floor(random))) = 4;

    
    Cells(1).NSC2Neuron = [];
    Cells(1).MC_division = [];
    
    
    %%% set the geometry
    geo(1).g = g;
    factor = params.factor;
    geo = MorphologyUpdate(geo,Cells,1,params,mi,factor);
    
    
    %%% find 3 bonds cells and increase their number of bonds to 4, if possible
    g1 = geo(1).g;
    for c = 1:length(g1.cells)-1
        if g1.dead(c) == 0
            bonds = g1.cells{c+1}; % the number of bonds of cell c 
            if bonds == 3
                [geo,StopSimulation] = AddBond(geo,Cells,num_of_initial_time_steps,c);
                if StopSimulation
                    return
                end
            end
        end
    end
    
    Initial_Cells = Cells;
    Initial_geo = geo;
end

