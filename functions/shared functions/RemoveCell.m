%% in case of 4 bonds
function [g1,new_state,StopSimulation] = RemoveCell(g,cell_num)
%REMOVECELL Will remove cell c from the lattice.
%   The function takes into account that cell c has only 4 bonds

    StopSimulation = 0;

    g1 = g;
    c = cell_num;                                               
    bonds = g1.cells{c+1};                                      % bonds = cells that has bonds with cell c
    verts_indx = g1.bonds(bonds,1);                             % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)
    verts_coordinates = getRelativePosition(g1,verts_indx,c);   % the coordinates of the verts
    
    if length(bonds) == 4
        
        neighbor1 = g1.bonds(bonds(1),4);
        neighbor2 = g1.bonds(bonds(2),4);
        neighbor3 = g1.bonds(bonds(3),4);
        neighbor4 = g1.bonds(bonds(4),4);
        neighbors_of_c = [ neighbor1 neighbor2 neighbor3 neighbor4 ];
        
        if neighbor1 == neighbor2 || neighbor1 == neighbor3 || neighbor1 == neighbor4 || neighbor2 == neighbor3 || neighbor2 == neighbor4 || neighbor3 == neighbor4
            disp('I need to figure out the cases and deal with them')
            keyboard
        else
            neighbors_of_neighbor1 = g1.bonds(g1.cells{neighbor1+1},4);
            neighbors_of_neighbor2 = g1.bonds(g1.cells{neighbor2+1},4);
            neighbors_of_neighbor3 = g1.bonds(g1.cells{neighbor3+1},4);
            neighbors_of_neighbor4 = g1.bonds(g1.cells{neighbor4+1},4);
            
            % find if any neighbor of cell c is also neighbor with adjacent neighbors of cell c 
            % (for example if neighbor1 is also neighbor with neighbor2 and neighbor4)
            neighbor1_with_neighbor2_and_neighbor4 = intersect(neighbors_of_neighbor1,[neighbor2 neighbor4]);
            neighbor2_with_neighbor1_and_neighbor3 = intersect(neighbors_of_neighbor2,[neighbor1 neighbor3]);
            neighbor3_with_neighbor2_and_neighbor4 = intersect(neighbors_of_neighbor3,[neighbor2 neighbor4]);
            neighbor4_with_neighbor1_and_neighbor3 = intersect(neighbors_of_neighbor4,[neighbor1 neighbor3]);
            
            %%% two options:
            % 1 - choose to delete 1 bond of cell c, in the condition that this bond's neighbor of c does not neighbor with the other neighbors of cell c
            % 2 - if option 1 is unavailable, delete 2 opposite bonds and connect between 2 opposite neighbors of c 
            
            % for option 1 - prefer the neighbor with the smallest area
            area1 = cellarea(g1,neighbor1);
            area2 = cellarea(g1,neighbor2);
            area3 = cellarea(g1,neighbor3);
            area4 = cellarea(g1,neighbor4);
            [~,AreaI] = sort([area1,area2,area3,area4]); % AreaI (index) is the order of neighbors of cell c from the smallest neighbor in area to the largest 
            
            NeighborL = NaN(1,4); % NeighborL (logical 0 or 1) is a logical vector for the 4 neighbors: 0 this neighbor of cell c has bond with other neighbors of cell c, 1 it has not 
            NeighborL(1) = isempty(neighbor1_with_neighbor2_and_neighbor4);
            NeighborL(2) = isempty(neighbor2_with_neighbor1_and_neighbor3);
            NeighborL(3) = isempty(neighbor3_with_neighbor2_and_neighbor4);
            NeighborL(4) = isempty(neighbor4_with_neighbor1_and_neighbor3);
                        
%             if isempty(neighbor1_with_neighbor2_and_neighbor4)
%                 delete_bond = bonds(1);
%                 neighbor = neighbor1;
            if NeighborL(AreaI(1))
                delete_bond = bonds(AreaI(1));
                neighbor = neighbors_of_c(AreaI(1));
            else
%                 if isempty(neighbor2_with_neighbor1_and_neighbor3)
%                     delete_bond = bonds(2);
%                     neighbor = neighbor2;
                if NeighborL(AreaI(2))
                    delete_bond = bonds(AreaI(2));
                    neighbor = neighbors_of_c(AreaI(2));
                else
%                     if isempty(neighbor3_with_neighbor2_and_neighbor4)
%                         delete_bond = bonds(3);
%                         neighbor = neighbor3;
                    if NeighborL(AreaI(3))
                        delete_bond = bonds(AreaI(3));
                        neighbor = neighbors_of_c(AreaI(3));
                    else
%                         if isempty(neighbor4_with_neighbor1_and_neighbor3)
%                             delete_bond = bonds(4);
%                             neighbor = neighbor4;
                        if NeighborL(AreaI(4))
                            delete_bond = bonds(AreaI(4));
                            neighbor = neighbors_of_c(AreaI(4));
                        else
                            delete_bond = 0;
                        end
                    end
                end
            end
            
            if delete_bond == 0  %% option 2
% %                 % turn cell c into NSC cell (cell state = 0)
% %                 new_state = 0;
                %%% delete two parallel bonds and join together two opposite neighbors 
                %%% (for example: delete bonds 2 and 4 and make neighbor 1 bond with neighbor 3)
                %%% choose to delete bonds from neighbors with maximum number of bonds
                number_of_bonds1 = length(neighbors_of_neighbor1);
                number_of_bonds2 = length(neighbors_of_neighbor2);
                number_of_bonds3 = length(neighbors_of_neighbor3);
                number_of_bonds4 = length(neighbors_of_neighbor4);
                
                if number_of_bonds1 < 4 && number_of_bonds2 < 4 || number_of_bonds2 < 4 && number_of_bonds3 < 4 || number_of_bonds3 < 4 && number_of_bonds4 < 4 || number_of_bonds4 < 4 && number_of_bonds1 < 4
                    disp('problem - two adjacent neighbors has less than 4 bonds')
                    if ~exist('new_state')
                        new_state = NaN;
                    end
                    StopSimulation = 1;
                    return
%                     keyboard
                end
                
                %%% if a neighbor has less than 4 bonds, choose not to reduce its number of bonds  
                less_than_4_bonds = find([number_of_bonds1 number_of_bonds2 number_of_bonds3 number_of_bonds4] < 4);
                if ~isempty(less_than_4_bonds)
                    if less_than_4_bonds(1) == 1 || less_than_4_bonds(1) == 3
                        bonds_to_remove = [2 4];    % the bonds to remove are bonds(bonds_to_remove)
                        bonds_to_keep = [1 3];      % the bonds to keep are bonds(bonds_to_keep)
                    else
                        bonds_to_remove = [1 3];
                        bonds_to_keep = [2 4];
                    end
                else
                    %%% if all neighbors have 4 or more bonds, choose to reduce the number of bonds from the neighbor with maximum number of bonds  
                    max_bonds = find([number_of_bonds1 number_of_bonds2 number_of_bonds3 number_of_bonds4] == max([number_of_bonds1 number_of_bonds2 number_of_bonds3 number_of_bonds4]));
                    if max_bonds(1) == 1 || max_bonds(1) == 3
                        bonds_to_remove = [1 3];
                        bonds_to_keep = [2 4];
                    else
                        bonds_to_remove = [2 4];
                        bonds_to_keep = [3 1];
                    end
                end
                % find the bonds number from the neighbors side, of the bonds to keep and remove   
                inverse_bonds_to_remove(1) = find(g1.bonds(:,3) == g1.bonds(bonds(bonds_to_remove(1)),4) & g1.bonds(:,4) == c);
                inverse_bonds_to_remove(2) = find(g1.bonds(:,3) == g1.bonds(bonds(bonds_to_remove(2)),4) & g1.bonds(:,4) == c);
                inverse_bonds_to_keep(1) = find(g1.bonds(:,3) == g1.bonds(bonds(bonds_to_keep(1)),4) & g1.bonds(:,4) == c);
                inverse_bonds_to_keep(2) = find(g1.bonds(:,3) == g1.bonds(bonds(bonds_to_keep(2)),4) & g1.bonds(:,4) == c);
                
                % set the new connection between the neighbors of c that are in "bonds_to_keep"
                % set neighbors
                g1.bonds(inverse_bonds_to_keep(1),4) = neighbors_of_c(bonds_to_keep(2)); 
                g1.bonds(inverse_bonds_to_keep(2),4) = neighbors_of_c(bonds_to_keep(1));
                % set vertices
                inverse_keep_bond_right = find(g1.bonds(:,2) == g1.bonds(inverse_bonds_to_keep(1),1) & g1.bonds(:,3) == neighbors_of_c(bonds_to_keep(1)));% set neighbors
                g1.bonds(inverse_keep_bond_right,2) = g1.bonds(inverse_bonds_to_keep(2),2);
                inverse_keep_bond_left = find(g1.bonds(:,1) == g1.bonds(inverse_bonds_to_keep(1),2) & g1.bonds(:,3) == neighbors_of_c(bonds_to_keep(1)));% set neighbors
                g1.bonds(inverse_keep_bond_left,1) = g1.bonds(inverse_bonds_to_keep(2),1);
                g1.bonds(inverse_bonds_to_keep(1),1) = g1.bonds(inverse_bonds_to_keep(2),2); 
                g1.bonds(inverse_bonds_to_keep(1),2) = g1.bonds(inverse_bonds_to_keep(2),1);
                
                % remove the bonds of "bonds_to_remove"
% %                 inverse_remove_bond_right = find(g1.bonds(:,2) == g1.bonds(inverse_bonds_to_remove(1),1) & g1.bonds(:,3) == neighbors_of_c(bonds_to_remove(1)));% set neighbors
% %                 g1.bonds(inverse_remove_bond_right,2) = g1.bonds(inverse_bonds_to_remove(1),2);
                vert_remove1 = g1.bonds(inverse_bonds_to_remove(1),1);
                vert_remain1 = g1.bonds(inverse_bonds_to_remove(1),2);
                inverse_remove_bonds11 = find(g1.bonds(:,1) == vert_remove1); % set neighbors
                inverse_remove_bonds12 = find(g1.bonds(:,2) == vert_remove1);
                g1.bonds(inverse_remove_bonds11,1) = vert_remain1;
                g1.bonds(inverse_remove_bonds12,2) = vert_remain1;
                g1.bonds(inverse_bonds_to_remove(1),:) = 0;
                neighbore_num1 = neighbors_of_c(bonds_to_remove(1));
                g1.cells{neighbore_num1+1} = g1.cells{neighbore_num1+1}(~ismember(g1.cells{neighbore_num1+1},inverse_bonds_to_remove(1)));
% %                 inverse_remove_bond_left = find(g1.bonds(:,1) == g1.bonds(inverse_bonds_to_remove(2),2) & g1.bonds(:,3) == neighbors_of_c(bonds_to_remove(2)));% set neighbors
% %                 g1.bonds(inverse_remove_bond_left,1) = g1.bonds(inverse_bonds_to_remove(2),1);
                vert_remove2 = g1.bonds(inverse_bonds_to_remove(2),2);
                vert_remain2 = g1.bonds(inverse_bonds_to_remove(2),1);
                inverse_remove_bonds21 = find(g1.bonds(:,1) == vert_remove2); % set neighbors
                inverse_remove_bonds22 = find(g1.bonds(:,2) == vert_remove2);
                g1.bonds(inverse_remove_bonds21,1) = vert_remain2;
                g1.bonds(inverse_remove_bonds22,2) = vert_remain2;
                g1.bonds(inverse_bonds_to_remove(2),:) = 0;
                neighbore_num2 = neighbors_of_c(bonds_to_remove(2));
                g1.cells{neighbore_num2+1} = g1.cells{neighbore_num2+1}(~ismember(g1.cells{neighbore_num2+1},inverse_bonds_to_remove(2)));
                
                % delete cell c 
                g1.bonds(bonds,:) = 0;
                g1.cells{c+1} = [];
                g1.dead(c) = 1;
                new_state = -1; % turn cell c into a dead/migrated cell (cell state = -1)
                
            else  %% option 1
                %%% delete only delete_bond
                % get the bonds numbers and verteces number
                delete_bond_verts = g1.bonds(delete_bond,1:2);
                adjacent_bond_right = bonds(g1.bonds(bonds,2) == delete_bond_verts(1));
                adjacent_bond_right_verts = g1.bonds(adjacent_bond_right,1:2);
                adjacent_bond_left = bonds(g1.bonds(bonds,1) == delete_bond_verts(2));
                adjacent_bond_left_verts = g1.bonds(adjacent_bond_left,1:2);
                opposite_bond = bonds(~ismember(bonds,[delete_bond,adjacent_bond_right,adjacent_bond_left]));
                opposite_bond_verts = g1.bonds(opposite_bond,1:2);
                
                inverse_bond = find(g1.bonds(:,1) == delete_bond_verts(2) & g1.bonds(:,2) == delete_bond_verts(1) & g1.bonds(:,4) == c);
                inverse_bond_right = find(g1.bonds(:,1) == adjacent_bond_right_verts(2) & g1.bonds(:,2) == adjacent_bond_right_verts(1) & g1.bonds(:,4) == c);
                inverse_bond_left = find(g1.bonds(:,1) == adjacent_bond_left_verts(2) & g1.bonds(:,2) == adjacent_bond_left_verts(1) & g1.bonds(:,4) == c);
                inverse_opposite_bond = find(g1.bonds(:,1) == opposite_bond_verts(2) & g1.bonds(:,2) == opposite_bond_verts(1) & g1.bonds(:,4) == c);
                
                % delete "delete_bond" and "inverse_bond"
                g1.bonds(delete_bond,:) = 0;
                g1.bonds(inverse_bond,:) = 0;
                
                % set the other bonds of cell c to be part of the bonds of cell "neighbor"
                neighbor_bonds = g1.cells{neighbor+1};
                inverse_bond_position = find(neighbor_bonds == inverse_bond);
                if inverse_bond_position == 1
                    g1.cells{neighbor+1} = [ adjacent_bond_left, opposite_bond, adjacent_bond_right, g1.cells{neighbor+1}(inverse_bond_position+1:end)];
                else
                    if inverse_bond_position == length(neighbor_bonds)
                        g1.cells{neighbor+1} = [g1.cells{neighbor+1}(1:inverse_bond_position-1), adjacent_bond_left, opposite_bond, adjacent_bond_right];
                    else
                        g1.cells{neighbor+1} = [g1.cells{neighbor+1}(1:inverse_bond_position-1), adjacent_bond_left, opposite_bond, adjacent_bond_right, g1.cells{neighbor+1}(inverse_bond_position+1:end)];
                    end
                end
                
                g1.bonds(adjacent_bond_right,3) = neighbor;
                g1.bonds(adjacent_bond_left,3) = neighbor;
                g1.bonds(opposite_bond,3) = neighbor;
                
                % update the inverse bonds
                g1.bonds(inverse_bond_right,4) = neighbor;
                g1.bonds(inverse_bond_left,4) = neighbor;
                g1.bonds(inverse_opposite_bond,4) = neighbor;
                
                % delete cell c 
                g1.cells{c+1} = [];
                g1.dead(c) = 1;
                new_state = -1; % turn cell c into a dead/migrated cell (cell state = -1)
                
            end
        end
    else
        disp('weird - num of bonds is not 4')
        if ~exist('new_state')
            new_state = NaN;
        end
        StopSimulation = 1;
        return
%         keyboard
    end
    
end

%% In case of 3 bonds
% % % function g1 = RemoveCell(g,cell_num)
% % % %REMOVECELL Will remove cell c from the lattice.
% % % %   The function takes into account that cell c has only 3 bonds
% % %     
% % %     g1 = g;
% % %     c = cell_num;                                               
% % %     bonds = g1.cells{c+1};                                      % bonds = cells that has bonds with cell_number
% % %     verts_indx = g1.bonds(bonds,1);                             % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)
% % %     verts_coordinates = getRelativePosition(g1,verts_indx,c);   % the coordinates of the verts
% % %     
% % %     if length(bonds) == 3
% % %         %%%% The idea is to delete 1 or 2 bonds of cell c, and eliminate this cell
% % %         %%%% there are 3 cases:
% % %         %%%% 
% % %         
% % %         neighbor1 = g1.bonds(bonds(1),4);
% % %         neighbor2 = g1.bonds(bonds(2),4);
% % %         neighbor3 = g1.bonds(bonds(3),4);
% % %         
% % %         %% case 1
% % %         if neighbor1 == neighbor2 && neighbor1 == neighbor3
% % %             disp('weird - cell c is within another cell')
% % %             keyboard
% % %         end
% % %         %% case 2
% % %         if neighbor1 == neighbor2 && neighbor1 ~= neighbor3
% % %             %%% delete the 2 bonds of cell c with neighbor1 (which is also neighbor2)
% % %             %%% update the third bond of cell c to be a new bond of neighbor1 
% % %             opposite_bonds = g1.cells{neighbor1+1}(g1.bonds(g1.cells{neighbor1+1},4) == c); % the bonds between neighbor1 and cell c
% % %             bond1_verts = g1.bonds(opposite_bonds(1),1:2);
% % %             bond2_verts = g1.bonds(opposite_bonds(2),1:2);
% % %             shared_vertex = intersect(bond1_verts,bond2_verts);
% % %             % update vertices of neighbor1 bond 1 to be the third bond of cell c (neighboring neighbor3) 
% % %             bond1_verts(bond1_verts == shared_vertex) = bond2_verts(bond2_verts ~= shared_vertex);
% % %             g1.bonds(opposite_bonds(1),1:2) = bond1_verts;
% % %             g1.bonds(opposite_bonds(1),4) = neighbor3;
% % %             % remove neighbor1 bond 2 from its bonds
% % %             g1.bonds(opposite_bonds(2),:) = 0;
% % %             g1.cells{neighbor1+1} = g1.cells{neighbor1 + 1}(~ismember(g1.cells{neighbor1+1},opposite_bonds(2)));
% % %             % update neighbor3 neighbor
% % %             bond3 = g1.cells{neighbor3+1}(g1.bonds(g1.cells{neighbor3+1},4) == c);
% % %             g1.bonds(bond3,4) = neighbor1;
% % %             % kill cell c and remove its bonds
% % %             g1.cells{c+1} = [];
% % %             g1.dead(c) = 1;
% % %             g1.bonds(bonds,:) = 0;
% % %             g1.verts(shared_vertex,:) = NaN;
% % %             % taking care of excess bond (if neighbor3 and neighbor1 have already shared a bond)
% % %             bonds_3to1 = g1.cells{neighbor3+1}(g1.bonds(g1.cells{neighbor3+1},4) == neighbor1);
% % %             if bonds_3to1 >= 2
% % %                 keyboard % check if it works properly
% % %                 g1 = ReduceBond(g1,bond3);
% % %             end
% % %         end
% % %         if neighbor1 == neighbor3 && neighbor1 ~= neighbor2
% % %             %%% delete the 2 bonds of cell c with neighbor1 (which is also neighbor3)
% % %             %%% update the third bond of cell c to be a new bond of neighbor1 
% % %             opposite_bonds = g1.cells{neighbor1+1}(g1.bonds(g1.cells{neighbor1+1},4) == c); % the bonds between neighbor1 and cell c
% % %             bond1_verts = g1.bonds(opposite_bonds(1),1:2);
% % %             bond2_verts = g1.bonds(opposite_bonds(2),1:2);
% % %             shared_vertex = intersect(bond1_verts,bond2_verts);
% % %             % update vertices of neighbor1 bond 1 to be the third bond of cell c (neighboring neighbor3) 
% % %             bond1_verts(bond1_verts == shared_vertex) = bond2_verts(bond2_verts ~= shared_vertex);
% % %             g1.bonds(opposite_bonds(1),1:2) = bond1_verts;
% % %             g1.bonds(opposite_bonds(1),4) = neighbor2;
% % %             % remove neighbor1 bond 2 from its bonds
% % %             g1.bonds(opposite_bonds(2),:) = 0;
% % %             g1.cells{neighbor1+1} = g1.cells{neighbor1 + 1}(~ismember(g1.cells{neighbor1+1},opposite_bonds(2)));
% % %             % update neighbor2 neighbor
% % %             bond3 = g1.cells{neighbor2+1}(g1.bonds(g1.cells{neighbor2+1},4) == c);
% % %             g1.bonds(bond3,4) = neighbor1;
% % %             % kill cell c and remove its bonds
% % %             g1.cells{c+1} = [];
% % %             g1.dead(c) = 1;
% % %             g1.bonds(bonds,:) = 0;
% % %             g1.verts(shared_vertex,:) = NaN;
% % %             % taking care of excess bond (if neighbor3 and neighbor1 have already shared a bond)
% % %             bonds_2to1 = g1.cells{neighbor2+1}(g1.bonds(g1.cells{neighbor2+1},4) == neighbor1);
% % %             if bonds_2to1 >= 2
% % %                 keyboard % check if it works properly
% % %                 g1 = ReduceBond(g1,bond3);
% % %             end
% % %         end
% % %         if neighbor2 == neighbor3 && neighbor2 ~= neighbor1
% % %             %%% delete the 2 bonds of cell c with neighbor2 (which is also neighbor3)
% % %             %%% update the third bond of cell c to be a new bond of neighbor2 
% % %             opposite_bonds = g1.cells{neighbor2+1}(g1.bonds(g1.cells{neighbor2+1},4) == c); % the bonds between neighbor1 and cell c
% % %             bond1_verts = g1.bonds(opposite_bonds(1),1:2);
% % %             bond2_verts = g1.bonds(opposite_bonds(2),1:2);
% % %             shared_vertex = intersect(bond1_verts,bond2_verts);
% % %             % update vertices of neighbor1 bond 1 to be the third bond of cell c (neighboring neighbor3) 
% % %             bond1_verts(bond1_verts == shared_vertex) = bond2_verts(bond2_verts ~= shared_vertex);
% % %             g1.bonds(opposite_bonds(1),1:2) = bond1_verts;
% % %             g1.bonds(opposite_bonds(1),4) = neighbor1;
% % %             % remove neighbor1 bond 2 from its bonds
% % %             g1.bonds(opposite_bonds(2),:) = 0;
% % %             g1.cells{neighbor2+1} = g1.cells{neighbor2 + 1}(~ismember(g1.cells{neighbor2+1},opposite_bonds(2)));
% % %             % update neighbor3 neighbor
% % %             bond3 = g1.cells{neighbor1+1}(g1.bonds(g1.cells{neighbor1+1},4) == c);
% % %             g1.bonds(bond3,4) = neighbor2;
% % %             % kill cell c and remove its bonds
% % %             g1.cells{c+1} = [];
% % %             g1.dead(c) = 1;
% % %             g1.bonds(bonds,:) = 0;
% % %             g1.verts(shared_vertex,:) = NaN;
% % %             % taking care of excess bond (if neighbor3 and neighbor1 have already shared a bond)
% % %             bonds_1to2 = g1.cells{neighbor1+1}(g1.bonds(g1.cells{neighbor1+1},4) == neighbor2);
% % %             if bonds_1to2 >= 2
% % %                 keyboard % check if it works properly
% % %                 g1 = ReduceBond(g1,bond3);
% % %             end
% % %         end
% % %         %% case 3
% % %         if neighbor1 ~= neighbor2 && neighbor1 ~= neighbor3 && neighbor2 ~= neighbor3
% % %             %%% delete the 1 bond of cell c with one neighbor upon these conditions: 
% % %             %%% The removal of the bond will not decrease the neighboring cell bonds from 3 to 2. 
% % %             %%% the prefared bond for removal will create, for its neighbor, bonds with new neigbors (1 or 2 new neighbors), 
% % %             %%% and, if possible, increase the number of bonds of this neighbor 
% % %             !!!!!!!!! incomplete
% % %         end
% % %     else
% % %         disp('weird - num of bonds is not 3')
% % %         keyboard
% % %     end
% % % end

