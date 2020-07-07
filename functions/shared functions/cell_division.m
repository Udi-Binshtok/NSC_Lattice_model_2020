function [g1, new_cell_num, weird] = cell_division( g, cell_num )
    % This code is taking into account that each cell has 3 or more bonds
    
    g1 = g;
    c = cell_num;
    c_new = length(g1.cells);
    bonds = g1.cells{c+1};                                      % bonds = cells that has bonds with cell_number
    verts_indx = g1.bonds(bonds,1);                             % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)
    verts_coordinates = getRelativePosition(g1,verts_indx,c);   % the coordinates of the verts as given by Micha + correction

    % if the cell's number of vertices is:
    % 6 or more - pick one vertex and one more at a distance of ~(#Of_Vertices/2) vertices away and divide the cell between these two
    % 5 - pick the longest bond and spilt it to 2, creating a new vertex. Pick another vertex at a distance of 3 vertices away from the new one and divide the cell between these two
    % 4 - pick the longest bond and also the opposing bond and split each of them to 2, creating 2 new vertices. Divide the cell between these two.
    % 3 - pick the longest bond and also the second longest bond and split each of them to 2, creating 2 new vertices. Divide the cell between these two.

    
    weird = 0;
    if length(bonds) <= 2
        weird = 1;          % this is for indicating an error
        new_cell_num = -1;  % this is for indicating an error as well
        return
    end
    
    %% first case (6 or more vertices)
    if length(bonds) >= 6
        vert1 = floor((length(bonds) - 1)*rand) + 1;    % choose one randon vertex
        vert2 = mod(vert1 + floor(length(verts_indx)/2),length(bonds));           % choose a second vertex at a distance of ~(#Of_Vertices/2) vertices from vert1
        if vert2 == 0
            vert2 = length(bonds);
        end
    end
    
    %% second case (5 vertices) - split one bond
    if length(bonds) == 5
        
        %%% calculating the length of each bond 
        for v = 1:(length(verts_indx) - 1) % v goes through the vertexes
            Delta_x = verts_coordinates(v,1) - verts_coordinates(v+1,1); % vertex v bonds with vertex v+1
            Delta_y = verts_coordinates(v,2) - verts_coordinates(v+1,2);
            bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        end
        v = length(verts_indx); % the last vertex bonds with the first one
        Delta_x = verts_coordinates(v,1) - verts_coordinates(1,1);
        Delta_y = verts_coordinates(v,2) - verts_coordinates(1,2);
        bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        
        %%% choose to split the longest bond into two halves (both the same size)   
        longest_bond_position = find(bonds_length == max(bonds_length));
        if length(longest_bond_position) > 1 % in the case their are several longest bonds at the same length
            longest_bond_position = longest_bond_position(1);
        end
        [g1,new_vert_number] = SplitBond(g1,cell_num,longest_bond_position);
        
        % get variabls again after change
        bonds = g1.cells{c+1};                                               
        verts_indx = g1.bonds(bonds,1);                                      
        verts_coordinates = getRelativePosition(g1,verts_indx,c);
        
        %%% choose vertices for cell division
        vert1 = find(verts_indx == new_vert_number);    % choose the first vertex to be the new one
        vert2 = mod(vert1 + 3,length(bonds));           % choose a second vertex at a distance of 3 vertices from vert1
        if vert2 == 0
            vert2 = length(bonds);
        end
    end
    
    %% third case (4 vertices) - split two bonds
    if length(bonds) == 4
        
        %%% 1st - split the longest bond
        %%% calculating the length of each bond 
        for v = 1:(length(verts_indx) - 1) % v goes through the vertexes
            Delta_x = verts_coordinates(v,1) - verts_coordinates(v+1,1); % vertex v bonds with vertex v+1
            Delta_y = verts_coordinates(v,2) - verts_coordinates(v+1,2);
            bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        end
        v = length(verts_indx); % the last vertex bonds with the first one
        Delta_x = verts_coordinates(v,1) - verts_coordinates(1,1);
        Delta_y = verts_coordinates(v,2) - verts_coordinates(1,2);
        bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        
        %%% choose to split the longest bond into two halves (both the same size)   
        longest_bond_position = find(bonds_length == max(bonds_length));
        if length(longest_bond_position) > 1 % in the case their are several longest bonds at the same length
            longest_bond_position = longest_bond_position(1);
        end
        [g1,new_vert1_number] = SplitBond(g1,cell_num,longest_bond_position);
        
        % get variabls again after change
        bonds = g1.cells{c+1};                                               
% % %         verts_indx = g1.bonds(bonds,1);                                      
% % %         verts_coordinates = getRelativePosition(g1,verts_indx,c);
        
        
% % %         %%% 2nd - split the 2nd longest bond
% % %         %%% calculating the length of each bond 
% % %         for v = 1:(length(verts_indx) - 1) % v goes through the vertexes
% % %             Delta_x = verts_coordinates(v,1) - verts_coordinates(v+1,1); % vertex v bonds with vertex v+1
% % %             Delta_y = verts_coordinates(v,2) - verts_coordinates(v+1,2);
% % %             bonds_lenght(v) = sqrt(Delta_x^2 + Delta_y^2);
% % %         end
% % %         v = length(verts_indx); % the last vertex bonds with the first one
% % %         Delta_x = verts_coordinates(v,1) - verts_coordinates(1,1);
% % %         Delta_y = verts_coordinates(v,2) - verts_coordinates(1,2);
% % %         bonds_lenght(v) = sqrt(Delta_x^2 + Delta_y^2);
% % %         
% % %         %%% choose to split the longest bond into two halves (both the same size)   
% % %         longest_bond_position = find(bonds_lenght == max(bonds_lenght));
% % %         % in the case where the longest bond is acctualy the two halves of the first bond that was spilt - I chhose to split the third longest bond
% % %         bonds_after_split1 = [ find(g1.bonds(bonds,1) == new_vert1_number) find(g1.bonds(bonds,2) == new_vert1_number) ];
% % % % %         if length(intersect(longest_bond_position,bonds_after_split1) == 1) == 2 
% % %         if length(intersect(longest_bond_position,bonds_after_split1)) == 2 
% % %             if longest_bond_position > 2 % this is the case where the longest bond is acctualy the two halves of the first bond that was spilt and they are the same length as other bonds
% % %                 bonds_lenght(bonds_after_split1) = bonds_lenght(bonds_after_split1) + [1,1];
% % %             end
% % %             [~,indx] = sort(bonds_lenght);
% % %             longest_bond_position = indx(3);
% % %         end
% % %         if length(longest_bond_position) > 1 % in the case their are several longest bonds at the same length
% % %             longest_bond_position = longest_bond_position(1);
% % %         end
% % %         [g1,new_vert2_number] = SplitBond(g1,cell_num,longest_bond_position);

        %%% 2nd - split the bond opposite of the longest bond (the bond that was split above). consider that now the cell has 5 bonds
        bonds_after_split1 = [ find(g1.bonds(bonds,2) == new_vert1_number) find(g1.bonds(bonds,1) == new_vert1_number) ];
        if length(intersect(bonds_after_split1,[1,5])) == 2
            opposit_to_longest_bond_position = 3;
        else
            if min(bonds_after_split1) == 1
                opposit_to_longest_bond_position = 4;
            else
                if min(bonds_after_split1) == 2
                    opposit_to_longest_bond_position = 5;
                else
                    opposit_to_longest_bond_position = min(bonds_after_split1) - 2;
                end
            end
        end        
        [g1,new_vert2_number] = SplitBond(g1,cell_num,opposit_to_longest_bond_position);
        
        % get variabls again after change
        bonds = g1.cells{c+1};                                               
        verts_indx = g1.bonds(bonds,1);                                      
% % %         verts_coordinates = getRelativePosition(g1,verts_indx,c);
        
        %%% choose vertices for cell division
        vert1 = find(verts_indx == new_vert2_number);   % choose the first vertex to be the new one
        vert2 = mod(vert1 + 3,length(bonds));           % choose a second vertex at a distance of 3 vertices from vert1
        if vert2 == 0
            vert2 = length(bonds);
        end
    end
    
    
    %% fourth case (3 vertices) - split two bonds
    if length(bonds) == 3
        
    % sanity check
    disp('please check if the cell division for 3 vertices is working properly')
%     keyboard
    
        %%% calculating the length of each bond 
        for v = 1:(length(verts_indx) - 1) % v goes through the vertexes
            Delta_x = verts_coordinates(v,1) - verts_coordinates(v+1,1); % vertex v bonds with vertex v+1
            Delta_y = verts_coordinates(v,2) - verts_coordinates(v+1,2);
            bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        end
        v = length(verts_indx); % the last vertex bonds with the first one
        Delta_x = verts_coordinates(v,1) - verts_coordinates(1,1);
        Delta_y = verts_coordinates(v,2) - verts_coordinates(1,2);
        bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
        
        %%% choose to split the longest bond into two halves (both the same size)  
        [~,ind] = sort(bonds_length);
        longest_bond_position = ind(end);
        second_longest_bond_position = ind(end-1);
        
        %%% 1st - split the longest bond
        [g1,new_vert1_number] = SplitBond(g1,cell_num,longest_bond_position);
                                                     
        %%% 2nd - split the second longest bond. consider that now the cell has 5 bonds
        [g1,new_vert2_number] = SplitBond(g1,cell_num,second_longest_bond_position);
        
        % get variabls again after change
        bonds = g1.cells{c+1};                                               
        verts_indx = g1.bonds(bonds,1);                                      
% % %         verts_coordinates = getRelativePosition(g1,verts_indx,c);
        
        %%% choose vertices for cell division
        vert1 = find(verts_indx == new_vert1_number);   % choose the first vertex to be the first new one
        vert2 = find(verts_indx == new_vert2_number);   % choose the second vertex to be the second new one
% % %         vert2 = mod(vert1 + 3,length(bonds));           % choose a second vertex at a distance of 3 vertices from vert1
        if vert2 == 0
            vert2 = length(bonds);
        end
    end
    
    
    %% cell division
    %%% add new cell
    start1 = find(g1.bonds(bonds,1) == verts_indx(vert1));
    end1 = find(g1.bonds(bonds,2) == verts_indx(vert2));
    if end1 < start1
        for i = 1:(length(bonds)-start1+1)
            bonds_cell_c(i) = bonds(start1 + i -1);
        end
        for i = 1:end1
            bonds_cell_c(length(bonds)-start1+1+i) = bonds(i);
        end
    else
        bonds_cell_c = bonds(start1:end1);
    end
    
    start2 = find(g1.bonds(bonds,1) == verts_indx(vert2));
    end2 = find(g1.bonds(bonds,2) == verts_indx(vert1));
    if end2 < start2
        for i = 1:(length(bonds)-start2+1)
            bonds_new_cell(i) = bonds(start2 + i -1);
        end
        for i = 1:end2
            bonds_new_cell(length(bonds)-start2+1+i) = bonds(i);
        end
    else
        bonds_new_cell = bonds(start2:end2);
    end
    
    % update g.cells with new cell
%     g1.cells{c+1} = g1.cells{c+1}(ismember(g1.cells{c+1},bonds_cell_c)); This is not always true (the one below is better right now)
    g1.cells{c+1} = bonds_cell_c;
    g1.cells{c_new + 1} = bonds_new_cell;
    
    %%% add new bond between cell c and the new cell
    Total_bonds_num = length(g1.bonds(:,1));
    % first in the direction of cell c to cell c_new
    g1.bonds(Total_bonds_num + 1,1) = verts_indx(vert2);
    g1.bonds(Total_bonds_num + 1,2) = verts_indx(vert1);
    g1.bonds(Total_bonds_num + 1,3) = c;
    g1.bonds(Total_bonds_num + 1,4) = c_new;
    g1.bonds(Total_bonds_num + 1,5) = 0.1;
    % second in the direction of cell c_new to cell c
    g1.bonds(Total_bonds_num + 2,1) = verts_indx(vert1);
    g1.bonds(Total_bonds_num + 2,2) = verts_indx(vert2);
    g1.bonds(Total_bonds_num + 2,3) = c_new;
    g1.bonds(Total_bonds_num + 2,4) = c;
    g1.bonds(Total_bonds_num + 2,5) = 0.1;
    
    % update g.cells with new bond
    d1 = Total_bonds_num + 1;
    d2 = Total_bonds_num + 2;
    g1.cells{c+1} = [ g1.cells{c+1} d1 ];                 
    g1.cells{c_new + 1} = [g1.cells{c_new+1} d2 ];       
    
    %%% update g.bonds - transfer bonds shared with cell c to cell c_new
    % direction of c_new to its neighbors
    g1.bonds(bonds_new_cell,3) = c_new;
    % direction of neighbors to c_new
    for b = 1:length(bonds_new_cell)     
        opposite_bond = find(g1.bonds(:,3) == g1.bonds(bonds_new_cell(b),4) & g1.bonds(:,4) == c); % bond number of the oppsite bond
        if length(opposite_bond) == 2 % it might be that the cell on opposite side of the bond shares two bonds with the original cell c
            % find which one of the opposite bonds is the opposite bond of c_new
            vertices_of_new_cell = g1.bonds(bonds_new_cell(b),1:2);
            vertices_of_opposing_cell1 = g1.bonds(opposite_bond(1),1:2);
            vertices_of_opposing_cell2 = g1.bonds(opposite_bond(2),1:2);
% % %             same_vertices(1) = length(intersect(vertices_of_new_cell,vertices_of_opposing_cell1) == 1);
% % %             same_vertices(2) = length(intersect(vertices_of_new_cell,vertices_of_opposing_cell2) == 1);
            same_vertices(1) = length(intersect(vertices_of_new_cell,vertices_of_opposing_cell1));
            same_vertices(2) = length(intersect(vertices_of_new_cell,vertices_of_opposing_cell2));
            opposite_bond = opposite_bond(find(same_vertices == 2));
            g1.bonds(opposite_bond,4) = c_new;
        else
            g1.bonds(opposite_bond,4) = c_new;
        end
    end
    
    %% Change the size of the new divided cells, according to some local avarage
    
    
% % %     % get the area of cell c and c_new after division 
% % %     divided_cells_area = [cellarea(g1,c), cellarea(g1,c_new)]; % areas of cells [c c_new] 
% % %     % get the area of the neighbors of cell c and c_new 
% % %     surrounding_cells = g.bonds(g.cells{c+1},4);
% % %     for i = 1:length(surrounding_cells)
% % %         neighbor_cells_area(i) = cellarea(g,surrounding_cells(i));
% % %     end 
% % %     
% % %     
% % %     % condition: if the average neighbor_cells_area is >= cell_division_increase_area_condition times the average divided_cells_area, then increase divded cells
% % %     cond = cell_division_increase_area_condition;
% % %     av_neighbor_cells_area = mean(neighbor_cells_area);
% % %     av_divided_cells_area = mean(divided_cells_area);
% % %     
% % %     % increase divided cells if needed
% % %     if av_neighbor_cells_area/av_divided_cells_area >= cond 
% % %         sd
% % %     end
% % %     
% % %     %%%%%%%% 
% % %     ratio_surr2div = [neighbor_cells_area./divided_cells_area(1) ; neighbor_cells_area./divided_cells_area(2)];
% % %     
% % %     
% % %     %%% choose the smaller cell to 
% % %     [~,indx] = sort(divided_cells_area);
% % %     if indx(1) == 1
% % %         enlarge_cell = c;
% % %     else
% % %         enlarge_cell = c_new;
% % %     end
% % %     vertices_cell_c = g.bonds(bonds_cell_c,1);
% % %     vertices_cell_c_new = g.bonds(bonds_new_cell,1);
% % %     
% % %     
% % %     %%% get the area of cell c before division
% % %     before_div_c_area = cellarea(g,c);
% % %     
    
    
    %% output
    % g1 = g1;
    new_cell_num = c_new;
    
    
    %% old, from micha
% % % %     nverts = size(g.verts,1);
% % % % 
% % % %     k = cell_num;
% % % % 
% % % %     nb = length(g.cells{k+1}); % number of bonds of cell k
% % % %     n1 = floor(rand*nb) + 1; % random num between 1 to number of bonds of cell k
% % % %     n2 = mod(n1 + floor(nb/2),nb);
% % % %     nn = sort([n1 n2]);
% % % % 
% % % %     % create divided cells
% % % %     c1 = [g.cells{k}(nn(2):nb) g.cells{k}(1:nn(1))];
% % % %     c2 = g.cells{k}(nn(1):nn(2));
% % % % 
% % % %     % old bonds
% % % %     b1 = g.cells{k+1}(nn(1));
% % % %     b2 = g.cells{k+1}(nn(2));
% % % % % %     b1 = g.cells{k}(nn(1));
% % % % % %     b2 = g.cells{k}(nn(2));
% % % % 
% % % %     % add new verts
% % % %     vx = 0.5*(g.verts(g.bonds(b1,1),1) + g.verts(g.bonds(b1,2),1));
% % % %     vy = 0.5*(g.verts(g.bonds(b1,1),2) + g.verts(g.bonds(b1,2),2));
% % % %     g.verts(nverts+1,:) = [vx vy 0];
% % % %     vx = 0.5*(g.verts(g.bonds(b2,1),1) + g.verts(g.bonds(b2,2),1));
% % % %     vy = 0.5*(g.verts(g.bonds(b2,1),2) + g.verts(g.bonds(b2,2),2));
% % % %     g.verts(nverts+2,:) = [vx vy 0];
% % % % 
% % % % 
% % % %     % update adjacent cells
% % % %     if g.bonds(b1,3) ~= k
% % % %         cn1 = g.bonds(b1,3);
% % % %     else
% % % %         cn1 = g.bonds(b1,4);
% % % %     end
% % % % 
% % % %     if g.bonds(b2,3) ~= k
% % % %         cn2 = g.bonds(b2,3);
% % % %     else
% % % %         cn2 = g.bonds(b2,4);
% % % %     end
end
    
