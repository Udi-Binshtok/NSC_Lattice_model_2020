function [geo1,StopSimulation] = ReduceVertexSharingCellNumber(geo,Cells,step,pert)
%ReduceVertexSharingCellNumber Will elimunate cases where cells are in touch with each
%other only by a vertex but with no bonds. The way it is done is by adding new
%bonds in places where a vertex has 4 or more cells that shares it, until there are only 3 cells that shares it.
%   pert = 0.05;
    
    StopSimulation = 0;
    
    geo1 = geo;
    g1 = geo1(step).g;
    Cells1 = Cells;
    
    aNSC = find(Cells1(step).states == 1); % which cells are aNSC at time step 'step'
    aNP = find(Cells1(step).states == 4); % which cells are aNP at time step 'step'
    div_aNP1 = find(Cells1(step).states == 4.1); % which aNP cells has already divided once and at state 4.1 at time step 'step'
    aNPs = [ aNP div_aNP1 ];
    
    vertices = g1.bonds(:,1);
    vertices = vertices(vertices ~= 0);
    vertices_index = unique(vertices); % these are all the vertices in time step 'step' 
    vertices_number_of_cells = histc(vertices, vertices_index); % this is the number of cells that shares the vertex 
    vertices = [ vertices_index , vertices_number_of_cells ];
    
    %% Find vertices that share 4 or more cells and reduce them to share 3 cells
    relevant_vertices = vertices(vertices_number_of_cells > 3,:);
    
    for v = 1:size(relevant_vertices,1)
        vertex_index = relevant_vertices(v,1);
        vertex_bonds = find(g1.bonds(:,1) == vertex_index); % which bonds have this vertex on their right edge 
        vertex_cells = g1.bonds(vertex_bonds,3); % which cells share this vertex 
        
        chosen_cells = NaN(1,size(vertex_cells,1) - 3); % the cells to be chosen for construction of new vertex and new bond  
        
        % choose aNP if possible
        shared_aNP = intersect(vertex_cells',aNPs); % find if any of the cells that share the vertex are aNP 
        shared_aNP = shared_aNP(randperm(length(shared_aNP)));
        if ~isempty(shared_aNP)
            c = 1;
            s = 1;
            while c <= size(chosen_cells,2) && s <= size(shared_aNP,2)
                chosen_cells(1,c) = shared_aNP(1,s);
                c = c + 1;
                s = s + 1;
            end
        end
        % choose aNSC if possible
        still_empty = find(isnan(chosen_cells(1,:)) == 1);
        if ~isempty(still_empty)
            shared_aNSC = intersect(vertex_cells',aNSC); % find if any of the cells that share the vertex are aNSC
            shared_aNSC = shared_aNSC(randperm(length(shared_aNSC)));
            if ~isempty(shared_aNSC)
                c = 1;
                s = 1;
                while c <= size(still_empty,2) && s <= size(shared_aNSC,2)
                    chosen_cells(1,still_empty(1,c)) = shared_aNSC(1,s);
                    c = c + 1;
                    s = s + 1;
                end
            end
        end
        % fill chosen cells with NSC if needed
        still_empty = find(isnan(chosen_cells(1,:)) == 1);
        if ~isempty(still_empty)
            shared_NSC = setdiff(setdiff(vertex_cells',aNPs),aNSC); % find all the cells that share the vertex that are not aNP and not aNSC
            shared_NSC = shared_NSC(randperm(length(shared_NSC)));
            c = 1;
            s = 1;
            while c <= size(still_empty,2) && s <= size(shared_NSC,2)
                chosen_cells(1,still_empty(1,c)) = shared_NSC(1,s);
                c = c + 1;
                s = s + 1;
            end
        end
        
        for i = 1:(size(vertex_cells,1) - 3) % reduce the number of shared cells to 3
            cell_number = chosen_cells(i);
            if isnan(cell_number)
                keyboard
            end
            [g1] = AddNewBond(g1,cell_number,vertex_index,pert);
        end
    end
    geo1(step).g = g1;
    
end

function [g1] = AddNewBond(g,cell_number,vertex_index,pert)
    
    g1 = g;

    %%% add a new vertex between the shared vertex and the center of the chosen cell that shares the vertex  
    cell_bonds = g.cells{cell_number+1}; % the bonds of the chosen cell                                              
    cell_verts_indx = g.bonds(cell_bonds,1); % the vertices of the chosen cell                                      
    cell_verts_coordinates = getRelativePosition(g,cell_verts_indx,cell_number);
    cell_CM = mean(cell_verts_coordinates); % the geometrical center of the chosen cell 
    vertex_coordinates_2 = g.verts(vertex_index,1:2);% the geometrical position of the shared vertex
    vertex_coordinates = cell_verts_coordinates(cell_verts_indx == vertex_index,:);% the geometrical position of the shared vertex
%     if vertex_coordinates_2 ~= vertex_coordinates
%         keyboard
%     end
    vector_shared_vert_to_cell = cell_CM - vertex_coordinates; % this is the vector from the shared vertex to chosen cell center
    
    % add new vertex
    new_vertex_number = size(g.verts,1) + 1;
    g1.verts(new_vertex_number,1:2) = vertex_coordinates + pert.*vector_shared_vert_to_cell; % the position of the new vertex is between the shared vertex and the center of the chosen cell
                                                                                          % pert will be the amount of movment along the shared_verrtex-cell_center axis
                                                                            
    % update the new vertex for the bonds of cell with its neighbors to the left and to the right of the shared vertex 
    cell_to_neighbor1_bond = cell_bonds(g.bonds(cell_bonds,2) == vertex_index); % neighbor 1 to the left of the shared vertex, relative to cell CM 
    cell_to_neighbor2_bond = cell_bonds(g.bonds(cell_bonds,1) == vertex_index); % neighbor 2 to the right of the shared vertex, relative to cell CM
    neighbor1_cellNumber = g.bonds(cell_to_neighbor1_bond,4);
    neighbor2_cellNumber = g.bonds(cell_to_neighbor2_bond,4);
    neighbor1_to_cell_bond = find(g.bonds(:,3) == neighbor1_cellNumber & g.bonds(:,4) == cell_number);
    neighbor2_to_cell_bond = find(g.bonds(:,3) == neighbor2_cellNumber & g.bonds(:,4) == cell_number);
    g1.bonds(cell_to_neighbor1_bond,2) = new_vertex_number; % update bond of cell with neighbor1 
    g1.bonds(cell_to_neighbor2_bond,1) = new_vertex_number; % update bond of cell with neighbor2
    g1.bonds(neighbor1_to_cell_bond,1) = new_vertex_number; % update bond of neighbor1 with cell 
    g1.bonds(neighbor2_to_cell_bond,2) = new_vertex_number; %#ok<FNDSB> % update bond of neighbor2 with cell
    
    % create a new bond between neighbor1 and neighbor2
    Total_bonds_num = size(g.bonds,1);
    % first in the direction of neighbor1 to neighbor2 
    neighbor1_to_neighbor2_bond = Total_bonds_num + 1;
    g1.bonds(neighbor1_to_neighbor2_bond,1) = vertex_index;
    g1.bonds(neighbor1_to_neighbor2_bond,2) = new_vertex_number;
    g1.bonds(neighbor1_to_neighbor2_bond,3) = neighbor1_cellNumber;
    g1.bonds(neighbor1_to_neighbor2_bond,4) = neighbor2_cellNumber;
    g1.bonds(neighbor1_to_neighbor2_bond,5) = 0.1;
    % next in the direction of neighbor2 to neighbor1 
    neighbor2_to_neighbor1_bond = Total_bonds_num + 2;
    g1.bonds(neighbor2_to_neighbor1_bond,1) = new_vertex_number;
    g1.bonds(neighbor2_to_neighbor1_bond,2) = vertex_index;
    g1.bonds(neighbor2_to_neighbor1_bond,3) = neighbor2_cellNumber;
    g1.bonds(neighbor2_to_neighbor1_bond,4) = neighbor1_cellNumber;
    g1.bonds(neighbor2_to_neighbor1_bond,5) = 0.1;

    % update g.cells with new bond
    % first the bonds of neighbor1
    neighbor1_to_cell_bondPosition = find(g.cells{neighbor1_cellNumber+1} == neighbor1_to_cell_bond); % find the position of bond from neighbor1 to cell within the bonds of neighbor1
    g1.cells{neighbor1_cellNumber+1} = [g.cells{neighbor1_cellNumber+1}(1:neighbor1_to_cell_bondPosition-1) neighbor1_to_neighbor2_bond g.cells{neighbor1_cellNumber+1}(neighbor1_to_cell_bondPosition:length(g.cells{neighbor1_cellNumber+1}))]; % add the bond from neighbor1 to neighbor 2 to the bonds of neighbor1, at the position of bond neighbor1 to cell
    % first the bonds of neighbor1
    neighbor2_to_cell_bondPosition = find(g.cells{neighbor2_cellNumber+1} == neighbor2_to_cell_bond); % find the position of bond from neighbor2 to cell within the bonds of neighbor2
    g1.cells{neighbor2_cellNumber+1} = [g.cells{neighbor2_cellNumber+1}(1:neighbor2_to_cell_bondPosition) neighbor2_to_neighbor1_bond g.cells{neighbor2_cellNumber+1}(neighbor2_to_cell_bondPosition+1:length(g.cells{neighbor2_cellNumber+1}))]; % add the bond from neighbor2 to neighbor 1 to the bonds of neighbor2, at the position of bond neighbor2 to cell
    
end

