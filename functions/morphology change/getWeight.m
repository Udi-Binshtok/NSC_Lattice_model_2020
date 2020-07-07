function we = getWeight(bo,g)
v1 = g.verts(g.bonds(bo,1),:);
v2 = g.verts(g.bonds(bo,2),:);
vec = v1-v2;
%we = 1+ abs(sin(atan2(vec(2),vec(1))));
we =1;
% if(g.bonds(bo,3)*g.bonds(bo,4)>1),
% %% if(g.dead(g.bonds(bo,3)) | g.dead(g.bonds(bo,4)))
% %%     we = 0.1;
% %%     return;
% %% end
%     if g.level(g.bonds(bo,3),1)*g.level(g.bonds(bo,4),1)<1,
%     if g.level(g.bonds(bo,3),1)+g.level(g.bonds(bo,4),1)>5,
%         we = 0.1;
%         %we  =1;
%     end
%     end
% 
% %    pos = 0.5*(g.vertex(g.bonds(bo,1))+g.vertex(g.bonds(bo,2)))
% %    we = pos(2)*pos(2)
% end
end