function dE = dHpot(g,i)
dE = zeros(2*length(g.verts),1);
vidx = g.bonds(g.cells{i+1},1);
dE(2*vidx) = (g.paras(3)*g.verts(vidx,2).^3)/length(vidx);

end