function [ cell_periodic_boundary_condition_location ] = Location_of_one_cell( g, cell_number )
%LOCATION function will give the location of a cell and its location due to periodic boundary condition
%   For a given g structure (by Micha Hirsh) this function will give the
%   coordinates of a cell in x and y and also polar, i.e 
%   cell_location=(x,y,r,teta)
%   
%   Also cell_periodic_boundary_condition_location is a tensor representing
%   the cell's location in the original lattice and the 8 "imaged" lattices surrounding the original lattice due to 2D
%   periodic boundary condition, meaning
%   cell_periodic_boundary_condition_location = 3x3x4 tensor, where 3x3 are
%   the lattices, for example (2,2) is the original lattice and (1,1) is
%   the upper-left imaged lattice. 4 is the position vector both kartizic and polar as described
%   above.


bonds = g.cells{cell_number + 1};                                   % bonds = cells that has bonds with cell_number
verts_indx = g.bonds(bonds,1);                                      % verts_indx = the inidicated number of one of the verts of each bond (all verts counted once that way)  
%verts_coordinates = g.verts(verts,1:2);                            % the coordinates of the verts as given by Micha
verts_coordinates = getRelativePosition(g,verts_indx,cell_number);  % the coordinates of the verts as given by Micha + correction
vector_average = sum(verts_coordinates,1)/length(verts_indx);       % the position of the center of the cell

x = vector_average(1);
y = vector_average(2);

R2 = (x)^2 + (y)^2;
R = sqrt(R2);

Teta = atan2(y,x);

cell_location = [x y R Teta];

% now periodic:

left = x - 2*pi;
right = x + 2*pi;
up = y + 2*pi;
down = y - 2*pi;

R11 = sqrt((left)^2+(up)^2);    Teta11 = atan(up/left);
R12 = sqrt((x)^2+(up)^2);       Teta12 = atan(up/x);
R13 = sqrt((right)^2+(up)^2);   Teta13 = atan(up/right);

R21 = sqrt((left)^2+(y)^2);     Teta21 = atan(y/left);
R23 = sqrt((right)^2+(y)^2);    Teta23 = atan(y/right);

R31 = sqrt((left)^2+(down)^2);  Teta31 = atan(down/left);
R32 = sqrt((x)^2+(down)^2);     Teta32 = atan(down/x);
R33 = sqrt((right)^2+(down)^2); Teta33 = atan(down/right);

cell_periodic_boundary_condition_location = NaN(3,3,4);

cell_periodic_boundary_condition_location(1,1,:) = [left up R11 Teta11];
cell_periodic_boundary_condition_location(1,2,:) = [x up R12 Teta12];
cell_periodic_boundary_condition_location(1,3,:) = [right up R13 Teta13];
cell_periodic_boundary_condition_location(2,1,:) = [left y R21 Teta21];
cell_periodic_boundary_condition_location(2,2,:) = cell_location;
cell_periodic_boundary_condition_location(2,3,:) = [right y R23 Teta23];
cell_periodic_boundary_condition_location(3,1,:) = [left down R31 Teta31];
cell_periodic_boundary_condition_location(3,2,:) = [x down R32 Teta32];
cell_periodic_boundary_condition_location(3,3,:) = [right down R33 Teta33];

end

