
function dE = dHlength(g,i)
    dE = zeros(2*length(g.verts),1);
    vidx = g.bonds(g.cells{i+1},1);
    vert = getRelativePosition(g,vidx,i);
    nb = length(vidx);
    for j = 1:nb
        prev = mod(j-2,nb)+1;
        next = mod(j,nb)+1;
        we = getWeight(g.cells{i+1}(nb),g);
        n1 = norm(vert(j,:)-vert(prev,:));
        n2 = norm(vert(j,:)-vert(next,:));
        if n1 > 0.0001
            dE(2*vidx(j)-1) = we*(vert(j,1) - vert(prev,1))/n1;
            dE(2*vidx(j)) = we* (vert(j,2) - vert(prev,2))/n1;
        else
            E(2*vidx(j)-1) = we;
            dE(2*vidx(j)) = we;
        end
 
        if n2 > 0.0001
            dE(2*vidx(j)-1) = dE(2*vidx(j)-1) + we*(vert(j,1)-vert(next,1))/n2;
            dE(2*vidx(j)) = dE(2*vidx(j)) + we* (vert(j,2)-vert(next,2))/n2;
        else
            E(2*vidx(j)-1) = dE(2*vidx(j)-1) + we;
            dE(2*vidx(j)) = dE(2*vidx(j)) + we;
        end
    end
    dE = dE* g.paras(2);
end
