function [g,StopSimulation] = kill_cells(g,cells)

StopSimulation = 0;

for ii = 1:length(cells)
    c = cells(ii);
    if g.dead(c) == 0
% %         while length(g.cells{c+1}) > 3
        while length(g.cells{c+1}) > 4
            vidx = g.bonds(g.cells{c+1},1);
            vert = getRelativePosition(g,vidx,c);
            nb = length(vidx);
            len = zeros(nb,1);
            for j = 1:nb
             len(j) = norm(vert(j,1:2)-vert(mod(j,nb)+1,1:2));
            end
            [ord, id] = sort(len);
% %             g = T1transition(g,g.cells{c+1}(id(1)));
            g = T1transition_Udi(g,g.cells{c+1}(id(1)),0.03);
        end
% %         g = T2transition(g,c);
        [g,~,StopSimulation] = RemoveCell(g,c);     % this function is instead of T2transition function
        if StopSimulation
            return
        end
    end
end