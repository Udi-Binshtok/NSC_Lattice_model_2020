function [g1,new_vert_number] = SplitBond(g,cell_num,bond_position)
%SPLITBOND will symmetrically split a bond of a cell into two halves (the bond of the opposing cell will also get split).
%   bond_position = out of the bonds of cell_num, in what position is the bond to be split (1sth bond? 2nd bond? ....)
    
    g1 = g;
    c = cell_num;
    bonds = g1.cells{c+1};                                     % bonds = cells that has bonds with cell_num
    verts_indx = g1.bonds(bonds,1);                            % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)
    verts_coordinates = getRelativePosition(g1,verts_indx,c);  % the coordinates of the verts as given by Micha + correction

    split.bond = bonds(bond_position);                              % bond number
    split.vert1_indx = verts_indx(bond_position);                   % first vert number
    split.vert1_coordinates = verts_coordinates(bond_position,:);   % first vert coordinate
    if bond_position == length(bonds) % the last vertex bonds with the first one
        split.vert2_indx = verts_indx(1);                                   % second vert number
        split.vert2_coordinates = verts_coordinates(1,:);                   % second vert coordinate
    else
        split.vert2_indx = verts_indx(bond_position + 1);
        split.vert2_coordinates = verts_coordinates(bond_position + 1,:);
    end
        
    %%% add new vertex (in the middle of the longest bond)
    new_vert_number = size(g1.verts,1) + 1;
    new_vert_coordinates = 0.5*(split.vert1_coordinates + split.vert2_coordinates);
    g1.verts(new_vert_number,:) = [new_vert_coordinates , 0]; 
        
    % update g.bonds - cut the bond into two halves
    %                  first half - from the 1st vertex to the new vertex (this will be updating the old bond)
    %                  second half - from the new vertex to the 2nd vertex (this will be a new bond)
        
    % old bond (first half)
    opposing_cell = g1.bonds(split.bond,4);
    opposing_bond = find(g1.bonds(:,3) == opposing_cell & g1.bonds(:,4) == c);
    g1.bonds(split.bond,2) = new_vert_number; % g1.bonds(split.bond,1) = split.vert1_indx;
    g1.bonds(opposing_bond,1) = new_vert_number; % g1.bonds(opposing_bond,2) = split.vert1_indx;
    % new bond (second half)
    Total_bonds_num = length(g1.bonds(:,1));
    % first in the direction of cell c to its neighbor
    g1.bonds(Total_bonds_num + 1,1) = new_vert_number;
    g1.bonds(Total_bonds_num + 1,2) = split.vert2_indx;
    g1.bonds(Total_bonds_num + 1,3) = c;
    g1.bonds(Total_bonds_num + 1,4) = opposing_cell;
    g1.bonds(Total_bonds_num + 1,5) = 0.1;
    % second in the direction of neighbor to cell c
    g1.bonds(Total_bonds_num + 2,1) = split.vert2_indx;
    g1.bonds(Total_bonds_num + 2,2) = new_vert_number;
    g1.bonds(Total_bonds_num + 2,3) = opposing_cell;
    g1.bonds(Total_bonds_num + 2,4) = c;
    g1.bonds(Total_bonds_num + 2,5) = 0.1;
        
    % update g.cells with new bond
    bonds(bond_position+1:end+1) = bonds(bond_position:end);
    bonds(bond_position+1) = Total_bonds_num + 1;
    g1.cells{c+1} = bonds;
    opposing_bonds = g1.cells{opposing_cell + 1};
    opposing_bond_position = find(opposing_bonds == opposing_bond);
    if opposing_bond_position == 1
        opposing_bonds(end+1) = Total_bonds_num + 2;
    else
        opposing_bonds(opposing_bond_position+1:end+1) = opposing_bonds(opposing_bond_position:end);
        opposing_bonds(opposing_bond_position) = Total_bonds_num + 2;
    end
    g1.cells{opposing_cell + 1} = opposing_bonds;
end

