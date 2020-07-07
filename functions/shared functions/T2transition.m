
%%% NOTE: there are two versions of T2transition here. You need to choose which of them to use
%%%       micha's version of T2transition is after my version.

% % % function g1 = T2transition(g,cell_num)
% % % %% NOTE: this is my version of T2transition (Micha's version is down)
% % %     
% % %     g1 = g;
% % %     c = cell_num;                                               
% % %     bonds = g1.cells{c+1};                                      % bonds = cells that has bonds with cell_number
% % %     verts_indx = g1.bonds(bonds,1);                             % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)
% % %     verts_coordinates = getRelativePosition(g1,verts_indx,c);   % the coordinates of the verts as given by Micha + correction
% % %     
% % %     if length(bonds) == 3
% % %         %%%% The idea is to split the longest bond of cell c into two halves and join the new vertex with the vertex that is opposite to it, at the center of cell c 
% % %         %%%% In that way eliminate cell c and keep the number of bonds of its neighbors the same as they are 
% % %         
% % %         %%% calculating the length of each bond 
% % %         for v = 1:(length(verts_indx) - 1)                                  % v goes through the vertexes
% % %             Delta_x = verts_coordinates(v,1) - verts_coordinates(v+1,1);    % vertex v bonds with vertex v+1
% % %             Delta_y = verts_coordinates(v,2) - verts_coordinates(v+1,2);
% % %             bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
% % %         end
% % %         v = length(verts_indx); % the last vertex bonds with the first one
% % %         Delta_x = verts_coordinates(v,1) - verts_coordinates(1,1);
% % %         Delta_y = verts_coordinates(v,2) - verts_coordinates(1,2);
% % %         bonds_length(v) = sqrt(Delta_x^2 + Delta_y^2);
% % %         
% % %         %%% choose to split the longest bond into two halves (both the same size)   
% % %         longest_bond_position = find(bonds_length == max(bonds_length));
% % %         if length(longest_bond_position) > 1 % in the case their are several longest bonds at the same length
% % %             longest_bond_position = longest_bond_position(1);
% % %         end
% % %         longest_bond_verts = g1.bonds(bonds(longest_bond_position),1:2);                % the vertices of the longest bond
% % %         opposite_vert_number = verts_indx(~ismember(verts_indx,longest_bond_verts));    % the opposite vertex to the longest bond
% % %         
% % %         [g1,new_vert_number] = SplitBond(g1,c,longest_bond_position);
% % %         
% % %         % get variabls again after change
% % %         bonds = g1.cells{c+1};                                               
% % %         verts_indx = g1.bonds(bonds,1);                                      
% % %         verts_coordinates = getRelativePosition(g1,verts_indx,c);
% % %         
% % %         
% % %         %%% join new vertex with opposite vertex at the center of cell c
% % %         % in this way the neighbor of cell c at the longest bond will become neighbor with the two cells at the other two bonds
% % %         CM = mean(verts_coordinates);               % geometrical center of cell c
% % %         g1.verts(opposite_vert_number,1:2) = CM;    % change the vertex coordinates of opposite_vert to be at the geometrical center of cell c
% % %         
% % %         % merge the halves of the splitted bond with the other bonds of cell c 
% % %         % first half bond:
% % %         bond_num = bonds(g1.bonds(bonds,1) == new_vert_number); % the bond num of one of the halves of the splitted bond
% % %         second_vert = g1.bonds(bond_num,2);                     % the second vertex of this half, the vertex that is not the new vertex 
% % %         neighbor = g1.bonds(bond_num,4);                        % cell c neighbor of this half bond
% % %         merge_bond = bonds(g1.bonds(bonds,1) == second_vert);   % cell c bond number of the bond that shares the second vertex
% % %         merge_neighbor = g1.bonds(merge_bond,4);                % cell c neighbor of the merge_bond
% % %         
% % %         opposite_neighbor_bond = find(g1.bonds(:,3) == neighbor & g1.bonds(:,1) == second_vert & g1.bonds(:,2) == new_vert_number);
% % %         g1.bonds(opposite_neighbor_bond,2) = opposite_vert_number;
% % %         g1.bonds(opposite_neighbor_bond,4) = merge_neighbor;
% % %         
% % %         opposite_merge_neighbor_bond = find(g1.bonds(:,3) == merge_neighbor & g1.bonds(:,1) == opposite_vert_number & g1.bonds(:,2) == second_vert);
% % %         g1.bonds(opposite_merge_neighbor_bond,4) = neighbor;
% % %         
% % %         g1.bonds(bond_num,:) = 0;
% % %         g1.bonds(merge_bond,:) = 0;
% % %         
% % %         % second half bond:
% % %         bond_num = bonds(g1.bonds(bonds,2) == new_vert_number);
% % %         second_vert = g1.bonds(bond_num,1);
% % %         neighbor = g1.bonds(bond_num,4);
% % %         merge_bond = bonds(g1.bonds(bonds,2) == second_vert);
% % %         merge_neighbor = g1.bonds(merge_bond,4);
% % %         
% % %         opposite_neighbor_bond = find(g1.bonds(:,3) == neighbor & g1.bonds(:,1) == new_vert_number & g1.bonds(:,2) == second_vert);
% % %         g1.bonds(opposite_neighbor_bond,1) = opposite_vert_number;
% % %         g1.bonds(opposite_neighbor_bond,4) = merge_neighbor;
% % %         
% % %         opposite_merge_neighbor_bond = find(g1.bonds(:,3) == merge_neighbor & g1.bonds(:,1) == second_vert & g1.bonds(:,2) == opposite_vert_number);
% % %         g1.bonds(opposite_merge_neighbor_bond,4) = neighbor;
% % %         
% % %         g1.bonds(bond_num,:) = 0;
% % %         g1.bonds(merge_bond,:) = 0;
% % %         
% % %         % delete cell c completely
% % %         g1.verts(new_vert_number,:) = []; % delete the new vertex (from the splitted bond of cell c)
% % %         g1.cells{c+1} = [];
% % %         g1.dead(c) = 1;
% % %     
% % %     else
% % %         disp('weird - num of bonds is not 3')
% % %         keyboard
% % %     end
       
function [ng, torem] = T2transition(g,ce)
%% This is micha's T2transition

    ng = g;
    if length(g.cells{ce+1}) <= 3               % if the number of bonds of cell ce is 3 or less
%         disp(['removing cell ' num2str(ce)]);
        % ng.cells{ce+1}
        bonds = ng.cells{ce+1};                 % the bonds of cell ce
        verts = ng.bonds(bonds,1);              % the indicated number of the vertices of cell ce
        pos = getRelativePosition(ng,verts,ce); % vertices' coordinates of cell ce (each row is a vertex, x column, y column)        
        v = mean(pos);                          % geometrical center of cell ce
        ng.verts(verts(1),1:2) = v;             % change the vertex coordinates of bond 1 of cell ce to be at the geometrical center of cell ce
        ng.bonds(ng.bonds(:,1) == verts(2),1) = verts(1);
        ng.bonds(ng.bonds(:,2) == verts(2),2) = verts(1);
        if length(bonds) == 3
            ng.bonds(ng.bonds(:,1) == verts(3),1) = verts(1);
            ng.bonds(ng.bonds(:,2) == verts(3),2) = verts(1);
        end
        for i = 1:length(bonds)
            bi = find(g.bonds(:,1) == g.bonds(bonds(i),2) & g.bonds(:,2) == g.bonds(bonds(i),1));
            bidx = 1;
            if length(bi) > 1 %% take the right one if many exist
                bidx = find(g.bonds(bi,4) == ce);
            end
            binv(i) = bi(bidx);
        end
        torem = [bonds binv];
        ng = remove_bond(ng,torem);
        ng.cells{ce+1} = [];
        ng.dead(ce) = 1;
    end
end