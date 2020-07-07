function g1 = ReduceBond(g,bond)
%REDUCEBOND removes bond "bond"(and its opposite bond) if the cell has 2 bonds with the same neighbor (bond is one of them) 
%   Detailed explanation goes here

    g1 = g;
    % find which of the bonds neighbor with the same cell 
    verts_indx = g1.bonds(bond,1:2);
    c = g1.bonds(bond,3);
    neighbor = g1.bonds(bond,4); 
    rigth_bond = find(g1.bonds(:,2) == verts_indx(1) & g1.bonds(:,3) == c);
    rigth_neighbor = g1.bonds(rigth_bond,4);
    left_bond = find(g1.bonds(:,1) == verts_indx(2) & g1.bonds(:,3) == c);
    left_neighbor = g1.bonds(left_bond,4);
    
    if rigth_neighbor == neighbor
        
        exceed_vert = verts_indx(1); 
        g1.bonds(rigth_bond,2) = verts_indx(2);
        g1.bonds(bond,:) = 0;
        g1.cells{c+1} = g.cells{c+1}(~ismember(g.cells{c+1},bond));
        
        opposite_bond = find(g1.bonds(:,1) == verts_indx(2) & g1.bonds(:,2) == exceed_vert & g1.bonds(:,3) == neighbor);
        g1.bonds(opposite_bond,:) = 0;
        opposite_rigth_bond = find(g1.bonds(:,1) == exceed_vert & g1.bonds(:,2) == g1.bonds(rigth_bond,1) & g1.bonds(:,3) == neighbor);
        g1.bonds(opposite_rigth_bond,1) = verts_indx(2);
        g1.cells{neighbor+1} = g.cells{neighbor+1}(~ismember(g.cells{neighbor+1},opposite_bond));
        
        g1.verts(exceed_vert,:) = NaN;
        
    else
        if left_neighbor == neighbor
            
            exceed_vert = verts_indx(2);
            g1.bonds(left_bond,1) = verts_indx(1);
            g1.bonds(bond,:) = 0;
            g1.cells{c+1} = g.cells{c+1}(~ismember(g.cells{c+1},bond));
            
            opposite_bond = find(g1.bonds(:,1) == exceed_vert & g1.bonds(:,2) == verts_indx(1) & g1.bonds(:,3) == neighbor);
            g1.bonds(opposite_bond,:) = 0;
            opposite_left_bond = find(g1.bonds(:,1) == g1.bonds(left_bond,2) & g1.bonds(:,2) == exceed_vert & g1.bonds(:,3) == neighbor);
            g1.bonds(opposite_left_bond,2) = verts_indx(1);
            g1.cells{neighbor+1} = g.cells{neighbor+1}(~ismember(g.cells{neighbor+1},opposite_bond));
            
            g1.verts(exceed_vert,:) = NaN;
            
        else
            disp('strange')
            keyboard
        end
    end
end

