function dE = denergy(v,gin)

    dE = zeros(2*length(gin.verts),1);
    g = insertverts(v,gin);
    for i = 1:length(g.cells)-1
        if g.dead(i) == 0
%          dEi = dHarea(g,i) + dHlength(g,i) + dHpot(g,i);
            dEi = dHarea(g,i) + dHlength_Roie(g,i);
            dE = dE + dEi;
        end
    end
    dE = dE';
end