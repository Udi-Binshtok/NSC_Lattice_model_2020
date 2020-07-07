function area = cellarea(g,cell_num)

% area=0;
c = cell_num;
verts_indx = g.bonds(g.cells{c+1},1);
verts_coordinates = getRelativePosition(g,verts_indx,c);
%%
% vert = vert*g.scale; 
area = polyarea(verts_coordinates(:,1),verts_coordinates(:,2));
