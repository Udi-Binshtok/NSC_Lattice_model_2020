function v = extractverts(g)
% extractverts arranges the vertices coordinates in a vector.
%   all the odd components of the vector are the x positions of the vertices 
%   all the even components of the vector are the y positions of the vertices 
    l = length(g.verts);
    v = zeros(2*l,1);
    v(1:2:2*l-1) = g.verts(:,1);
    v(2:2:2*l) = g.verts(:,2);
% for i=1:length(g.verts)
%     v(2*i-1:2*i)=g.verts(i,1:2)';
% end