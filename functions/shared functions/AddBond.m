function [geo1,StopSimulation] = AddBond(geo,Cells,step,cell_num)
%ADDBOND Will add a bond to a cell (can be by using T1transition)
%   The added bond to cell 'cell_num' can reduce or add a bond from a different cell.
%   [The idea is to reduce a bond from a cell which has a large number of
%   bonds, and more then 4 bonds, and by that increase the number of bonds
%   of 'cell_num' (which, for example, might have 3 bonds at the moment)]
    
    StopSimulation = 0;
    
    c = cell_num;
    geo1 = geo;
    g1 = geo1(step).g;
    Cells1 = Cells;
    
    bonds = g1.cells{c+1};                          % bonds = the bonds of cell c
    verts_indx = g1.bonds(bonds,1);                 % verts_indx = the inidicated number of one of the vertices of each of the bonds of cell c (all the verts of cell c are counted once that way)
%   vidx = g.bonds(g.cells{c+1},1);
%     vert = getRelativePosition(g1,verts_indx,c);    % vert = vertices' coordinates of cell c (each row is a vertex, x column, y column)
    
    
    %%% Getting the number of neighbors sharing each vertex of cell c
    verts_sharing = NaN(length(verts_indx),1);
    for i = 1:length(verts_indx)
        verts_sharing(i) = length(find(g1.bonds(:,1) == verts_indx(i))); % the number of cells sharing each of the vertices of cell c
    end
    possible_vertices = find(verts_sharing >= 4); % find vertices of cell c that shares at least 4 neighbors
    
    %%% Two options for adding a bond to cell c:
    if ~isempty(possible_vertices)
        
    %% option 1: 4 (or more) cells shares 1 vertex - add new bond (in principle the lattice should not have these cases of more than 3 cells per vertex) 
        %%% Choose to add a bond for cell c with an opposite cell, cell c_opposite, that shares a vertex in "possible_vertex" 
        %%% Cell c_opposite is not a neighbor of cell c. Try to avoid choosing c_opposite that is a Neuron 
        %%% In that way both cell c and cell c_opposite increases the number of their bonds 
        neighbors_of_c = g1.bonds(bonds,4);
        if length(neighbors_of_c) > 3
            disp('at the moment, the function AddBond is working on cells with 3 bonds only')
            StopSimulation = 1;
            return
%             keyboard
        end
        Neurons = find(Cells1(step).states == 5);     % which cells are Neurons at time step 'step'
        for v = 1:length(possible_vertices)
            cells_sharing_vertex = g1.bonds(g1.bonds(:,1) == verts_indx(possible_vertices(v)),3);
            cells_sharing_vertex_without_c = cells_sharing_vertex(~ismember(cells_sharing_vertex,c));
            opposite_cells = cells_sharing_vertex_without_c(~ismember(cells_sharing_vertex_without_c,neighbors_of_c)); % these are cells that shares the vertex v in "possible_vertices" and that are not neighbors of c 
            
            c_opposite = opposite_cells(1); % this is a choise for cell c_opposite, but it can still be a Neuron
            shared_vert = verts_indx(possible_vertices(v)); % this is the shared vertex between cell c and cell c_opposite
            
            opposite_cells_not_neurons = opposite_cells(~ismember(opposite_cells,Neurons)); % these are cells that shares the vertex v in "possible_vertices" , that are not neighbors of c and that are not neurons
            if ~isempty(opposite_cells_not_neurons)
                c_opposite = opposite_cells_not_neurons(1); % this is a choise for cell c_opposite, that is also not a Neuron (break the loop)
                shared_vert = verts_indx(possible_vertices(v)); % (break the loop)
                break
            end
        end
        
        bonds_of_c_opposite = g1.cells{c_opposite+1}; % the bonds of cell c_opposite
        neighbors_of_c_opposite = g1.bonds(bonds_of_c_opposite,4);
        shared_neighbors = intersect(neighbors_of_c,neighbors_of_c_opposite); % these are neighbors of both c and c_opposite
        
        %%% two possibilities - cell c and cell c_opposite shares a neighbor or not
        if ~isempty(shared_neighbors)
            %%% if cell c_opposite shares a neighbor with cell c, add new vertex between shared_vert and the center of the shared neighbor 
            neighbor = shared_neighbors(1); % choose one of the shared neighbors in case there is more than one shared neighbor
            neighbor_bonds = g1.cells{neighbor+1}; % the bonds of the shared neighbor                                              
            neighbor_verts_indx = g1.bonds(neighbor_bonds,1);                                      
            neighbor_verts_coordinates = getRelativePosition(g1,neighbor_verts_indx,neighbor);
            neighbor_CM = mean(neighbor_verts_coordinates); % the geometrical center of the shared neighbor 
            shared_vert_coordinates = g1.verts(shared_vert,1:2);% the geometrical position of the shared vertex (between cell c, cell c_opposite and their shared neighbor) 
            vector_shared_vert_to_neighbor = neighbor_CM - shared_vert_coordinates; % this is the vector from shared_vert to neighbor center
            
            % add new vertex
            new_vert_number = size(g1.verts,1) + 1;
            g1.verts(new_vert_number,1:2) = shared_vert_coordinates + 0.05*vector_shared_vert_to_neighbor; % the position of the new vertex is between the shared vertex and the center of the shared neighbor
            
            % rearrange bonds between c, c_opposite and shared neighbor
            % (take their shared vertex to be the new vertex instead )
            bond_c2neighbor = find(g1.bonds(:,3) == c & g1.bonds(:,4) == neighbor);
            shared_vert_pos_c2neighbor = find(g1.bonds(bond_c2neighbor,1:2) == shared_vert);
            g1.bonds(bond_c2neighbor,shared_vert_pos_c2neighbor) = new_vert_number; % set the shared vertex of bond c2neighbor to be the new vertex

            bond_neighbor2c = find(g1.bonds(:,3) == neighbor & g1.bonds(:,4) == c);
            shared_vert_pos_neighbor2c = find(g1.bonds(bond_neighbor2c,1:2) == shared_vert);
            g1.bonds(bond_neighbor2c,shared_vert_pos_neighbor2c) = new_vert_number; % set the shared vertex of bond neighbor2c to be the new vertex
            
            bond_c_opposite2neighbor = find(g1.bonds(:,3) == c_opposite & g1.bonds(:,4) == neighbor);
            shared_vert_pos_c_opposite2neighbor = find(g1.bonds(bond_c_opposite2neighbor,1:2) == shared_vert);
            g1.bonds(bond_c_opposite2neighbor,shared_vert_pos_c_opposite2neighbor) = new_vert_number; % set the shared vertex of bond c_opposite2neighbor to be the new vertex
            
            bond_neighbor2c_opposite = find(g1.bonds(:,3) == neighbor & g1.bonds(:,4) == c_opposite);
            shared_vert_pos_neighbor2c_opposite = find(g1.bonds(bond_neighbor2c_opposite,1:2) == shared_vert);
            g1.bonds(bond_neighbor2c_opposite,shared_vert_pos_neighbor2c_opposite) = new_vert_number; % set the shared vertex of bond neighbor2c_opposite to be the new vertex
            
            % new bond
            Total_bonds_num = length(g1.bonds(:,1));
            % first in the direction of cell c to c_opposite
            if shared_vert_pos_c2neighbor == 1
                second_vert_pos_c2neighbor = 2;
            else
                second_vert_pos_c2neighbor = 1;
            end
            g1.bonds(Total_bonds_num + 1,shared_vert_pos_c2neighbor) = shared_vert;
            g1.bonds(Total_bonds_num + 1,second_vert_pos_c2neighbor) = new_vert_number;
            g1.bonds(Total_bonds_num + 1,3) = c;
            g1.bonds(Total_bonds_num + 1,4) = c_opposite;
            g1.bonds(Total_bonds_num + 1,5) = 0.1;
            % second in the direction of c_opposite to cell c
            if shared_vert_pos_c_opposite2neighbor == 1
                second_vert_pos_c_opposite2neighbor = 2;
            else
                second_vert_pos_c_opposite2neighbor = 1;
            end
            g1.bonds(Total_bonds_num + 2,shared_vert_pos_c_opposite2neighbor) = shared_vert;
            g1.bonds(Total_bonds_num + 2,second_vert_pos_c_opposite2neighbor) = new_vert_number;
            g1.bonds(Total_bonds_num + 2,3) = c_opposite;
            g1.bonds(Total_bonds_num + 2,4) = c;
            g1.bonds(Total_bonds_num + 2,5) = 0.1;
            
            % update g.cells with new bond
            % first the bonds of cell c
            c_bond_position = find(bonds == bond_c2neighbor);
            if shared_vert_pos_c2neighbor == 1
                plus = 0;
            else
                plus = 1;
            end
            new_bonds_c = bonds;
            if c_bond_position == length(new_bonds_c) && plus == 1
                new_bonds_c(end+plus) = Total_bonds_num + 1;
            else
                new_bonds_c(c_bond_position+1+plus:end+1) = new_bonds_c(c_bond_position+plus:end);
                new_bonds_c(c_bond_position+plus) = Total_bonds_num + 1;
            end
            g1.cells{c+1} = new_bonds_c;
            % second the bonds of cell c_opposite
            c_opposite_bonds = g1.cells{c_opposite + 1};
            c_opposite_bond_position = find(c_opposite_bonds == bond_c_opposite2neighbor);
            if shared_vert_pos_c_opposite2neighbor == 1
                plus = 0;
            else
                plus = 1;
            end
            new_bonds_c_opposite = c_opposite_bonds;
            if c_opposite_bond_position == length(new_bonds_c_opposite) && plus == 1
                new_bonds_c_opposite(end+plus) = Total_bonds_num + 2;
            else
                new_bonds_c_opposite(c_opposite_bond_position+1+plus:end+1) = new_bonds_c_opposite(c_opposite_bond_position+plus:end);
                new_bonds_c_opposite(c_opposite_bond_position+plus) = Total_bonds_num + 2;
            end
            g1.cells{c_opposite+1} = new_bonds_c_opposite;
            
        else
            %%% if cell c_opposite does not share a neighbor with cell c, add new vertex between shared_vert and the center of the shared neighbor 
            % first - map and count the cells between cell c and cell c_opposite, from both sides 
            % Clockwise - first while loop from cell c to cell c_opposite 
            %             second while loop from cell c_opposite to cell c 
            current_bond = bonds(g1.bonds(bonds,1) == shared_vert);
            cells_c2c_opposite = [];
            flag1 = 1;
            while flag1
                n = g1.bonds(current_bond,4); % cell number of the neighbor with current_bond
                if n == c_opposite
                    flag1 = 0;
                    break
                end
                cells_c2c_opposite = [cells_c2c_opposite n];
                bonds_of_n = g1.cells{n+1};
                current_bond = bonds_of_n(g1.bonds(bonds_of_n,1) == shared_vert);
            end
            bonds_of_c_opposite = g1.cells{c_opposite+1};
            current_bond = bonds_of_c_opposite(g1.bonds(bonds_of_c_opposite,1) == shared_vert);
            cells_c_opposite2c = [];
            flag2 = 1;
            while flag2
                n = g1.bonds(current_bond,4); % cell number of the neighbor with current_bond
                if n == c
                    flag2 = 0;
                    break
                end
                cells_c_opposite2c = [cells_c_opposite2c n];
                bonds_of_n = g1.cells{n+1};
                current_bond = bonds_of_n(g1.bonds(bonds_of_n,1) == shared_vert);
            end
            % new vertex at the side of the least cells between cell c and c_opposite
            % the new vertex will be between the shared vertex position and the center of the neighbor of cell c at that side 
            if length(cells_c2c_opposite) < length(cells_c_opposite2c)
                n_cells = cells_c2c_opposite;
                neighbor = cells_c2c_opposite(1);
                c_opposite_neighbor = cells_c2c_opposite(end);
% %                 shared_vert_pos_c2neighbor = 1;
% %                 shared_vert_pos_c_opposite2neighbor = 2;
            else
                n_cells = cells_c_opposite2c;
                neighbor = cells_c_opposite2c(end);
                c_opposite_neighbor = cells_c_opposite2c(1);
% %                 shared_vert_pos_c2neighbor = 2;
% %                 shared_vert_pos_c_opposite2neighbor = 1;
            end
            
            neighbor_bonds = g1.cells{neighbor+1}; % the bonds of the chosen neighbor                                              
            neighbor_verts_indx = g1.bonds(neighbor_bonds,1);                                      
            neighbor_verts_coordinates = getRelativePosition(g1,neighbor_verts_indx,neighbor);
            neighbor_CM = mean(neighbor_verts_coordinates); % the geometrical center of the chosen neighbor 
            shared_vert_coordinates = g1.verts(shared_vert,1:2);% the geometrical position of the shared vertex (between cell c, cell c_opposite and the chosen neighbor) 
            vector_shared_vert_to_neighbor = neighbor_CM - shared_vert_coordinates; % this is the vector from shared_vert to neighbor center
            
            % add new vertex
            new_vert_number = size(g1.verts,1) + 1;
            g1.verts(new_vert_number,1:2) = shared_vert_coordinates + 0.05*vector_shared_vert_to_neighbor; % the position of the new vertex is between the shared vertex and the center of the chosen neighbor
            
            % rearrange bonds between c, c_opposite and shared neighbor
            % (take their shared vertex to be the new vertex instead )
            bond_c2neighbor = find(g1.bonds(:,3) == c & g1.bonds(:,4) == neighbor);
            shared_vert_pos_c2neighbor = find(g1.bonds(bond_c2neighbor,1:2) == shared_vert);
            g1.bonds(bond_c2neighbor,shared_vert_pos_c2neighbor) = new_vert_number; % set the shared vertex of bond c2neighbor to be the new vertex

% %             bond_neighbor2c = find(g1.bonds(:,3) == neighbor & g1.bonds(:,4) == c);
% %             shared_vert_pos_neighbor2c = find(g1.bonds(bond_neighbor2c,1:2) == shared_vert);
% %             g1.bonds(bond_neighbor2c,shared_vert_pos_neighbor2c) = new_vert_number; % set the shared vertex of bond neighbor2c to be the new vertex
            
            bond_c_opposite2neighbor = find(g1.bonds(:,3) == c_opposite & g1.bonds(:,4) == c_opposite_neighbor);
            shared_vert_pos_c_opposite2neighbor = find(g1.bonds(bond_c_opposite2neighbor,1:2) == shared_vert);
            g1.bonds(bond_c_opposite2neighbor,shared_vert_pos_c_opposite2neighbor) = new_vert_number; % set the shared vertex of bond c_opposite2neighbor to be the new vertex
            
% %             bond_neighbor2c_opposite = find(g1.bonds(:,3) == neighbor & g1.bonds(:,4) == c_opposite);
% %             shared_vert_pos_neighbor2c_opposite = find(g1.bonds(bond_neighbor2c_opposite,1:2) == shared_vert);
% %             g1.bonds(bond_neighbor2c_opposite,shared_vert_pos_neighbor2c_opposite) = new_vert_number; % set the shared vertex of bond neighbor2c_opposite to be the new vertex
            
            for i = 1:length(n_cells)
                n = n_cells(i);
                n_bonds = g1.cells{n+1};
                n_relevant_bond1 = n_bonds(g1.bonds(n_bonds,1) == shared_vert);
                n_relevant_bond2 = n_bonds(g1.bonds(n_bonds,2) == shared_vert);
                g1.bonds(n_relevant_bond1,1) = new_vert_number; % set the shared vertex of n_relevant_bond1 to be the new vertex
                g1.bonds(n_relevant_bond2,2) = new_vert_number; % set the shared vertex of n_relevant_bond2 to be the new vertex
            end
            
            % new bond
            Total_bonds_num = length(g1.bonds(:,1));
            % first in the direction of cell c to c_opposite
            if shared_vert_pos_c2neighbor == 1
                second_vert_pos_c2neighbor = 2;
            else
                second_vert_pos_c2neighbor = 1;
            end
            g1.bonds(Total_bonds_num + 1,shared_vert_pos_c2neighbor) = shared_vert;
            g1.bonds(Total_bonds_num + 1,second_vert_pos_c2neighbor) = new_vert_number;
            g1.bonds(Total_bonds_num + 1,3) = c;
            g1.bonds(Total_bonds_num + 1,4) = c_opposite;
            g1.bonds(Total_bonds_num + 1,5) = 0.1;
            % second in the direction of c_opposite to cell c
            if shared_vert_pos_c_opposite2neighbor == 1
                second_vert_pos_c_opposite2neighbor = 2;
            else
                second_vert_pos_c_opposite2neighbor = 1;
            end
            g1.bonds(Total_bonds_num + 2,shared_vert_pos_c_opposite2neighbor) = shared_vert;
            g1.bonds(Total_bonds_num + 2,second_vert_pos_c_opposite2neighbor) = new_vert_number;
            g1.bonds(Total_bonds_num + 2,3) = c_opposite;
            g1.bonds(Total_bonds_num + 2,4) = c;
            g1.bonds(Total_bonds_num + 2,5) = 0.1;
            
            % update g.cells with new bond
            % first the bonds of cell c
            c_bond_position = find(bonds == bond_c2neighbor);
            if shared_vert_pos_c2neighbor == 1
                plus = 0;
            else
                plus = 1;
            end
            new_bonds_c = bonds;
            if c_bond_position == length(new_bonds_c) && plus == 1
                new_bonds_c(end+plus) = Total_bonds_num + 1;
            else
                new_bonds_c(c_bond_position+1+plus:end+1) = new_bonds_c(c_bond_position+plus:end);
                new_bonds_c(c_bond_position+plus) = Total_bonds_num + 1;
            end
            g1.cells{c+1} = new_bonds_c;
            % second the bonds of cell c_opposite
            c_opposite_bonds = g1.cells{c_opposite + 1};
            c_opposite_bond_position = find(c_opposite_bonds == bond_c_opposite2neighbor);
            if shared_vert_pos_c_opposite2neighbor == 1
                plus = 0;
            else
                plus = 1;
            end
            new_bonds_c_opposite = c_opposite_bonds;
            if c_opposite_bond_position == length(new_bonds_c_opposite) && plus == 1
                new_bonds_c_opposite(end+plus) = Total_bonds_num + 2;
            else
                new_bonds_c_opposite(c_opposite_bond_position+1+plus:end+1) = new_bonds_c_opposite(c_opposite_bond_position+plus:end);
                new_bonds_c_opposite(c_opposite_bond_position+plus) = Total_bonds_num + 2;
            end
            g1.cells{c_opposite+1} = new_bonds_c_opposite;
            
        
        end
        
         
        
    else
    %% option 2: 3 cells shares 1 vertex - apply T1transition
        % In T1transition cell c will gain a bond but two of its neighbors will loose a bond,
        % therefore the neighbors that looses a bond should have 5 or more bonds so they will end up with 4 or more bonds 
        neighbors_of_c = g1.bonds(bonds,4);
        if length(neighbors_of_c) > 3
            disp('at the moment, the function AddBond is working on cells with 3 bonds only')
            StopSimulation = 1;
            return
%             keyboard
        end
        num_of_bonds = NaN(length(neighbors_of_c),1);
        for i = 1:length(neighbors_of_c)
            num_of_bonds(i) = length(g1.cells{neighbors_of_c(i)+1});
        end
        possible_neighbors = neighbors_of_c(num_of_bonds >= 5);
        if length(possible_neighbors) <= 1
            disp('problem - sorrounding cells of cell c have 4 or less bonds')
            StopSimulation = 1;
            return
%             keyboard
        else
            possible_neighbors_num_of_bonds = num_of_bonds(num_of_bonds >= 5);
            [~,ind] = sort(possible_neighbors_num_of_bonds,'descend');
            chosen_neighbors = possible_neighbors(ind(1:2));    
        end
        bonds_of_chosen_neighbor1 = g1.cells{chosen_neighbors(1)+1};
        chosen_bond = bonds_of_chosen_neighbor1(g1.bonds(bonds_of_chosen_neighbor1,4) == chosen_neighbors(2));
        g1 = T1transition_Udi(g1,chosen_bond,0.03);
    end
    
    geo1(step).g = g1;
end

