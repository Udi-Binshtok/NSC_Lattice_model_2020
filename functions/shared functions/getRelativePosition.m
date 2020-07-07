function pos = getRelativePosition(g,verts_indx,cell_number)
% v are the vertexes of cell i.
pos = g.verts(verts_indx,1:2);  % vertices coordinates of cell_number. pos is (#_Of_Verts)X2 matrix (x column, y column)
if(g.bc == 1)                   % bc = 1 means periodic boundary condition
    for d = 1:2                 % x and y coordinates 
        ap = g.verts(g.bonds(g.cells{cell_number+1},1),d);  % all vertices (positions) of cell_number
        [~, pid] = min(abs(ap));                            % p = value ; pid = posision in vector ap (it finds the closest vert to zero in ap)
        p = ap(pid);                                        % p = the vert closest to zero in ap
%         s = g.scale(1,1);
        s = 1;
        mup = find(abs(pos(:,d)+s*2*pi-p) < abs(pos(:,d)-p));
        mdown = find(abs(pos(:,d)-s*2*pi-p) < abs(pos(:,d)-p));
        pos(mup,d) = pos(mup,d)+(s*2*pi);                     % This set the vertexes coordinates of cell_number in the correct position, due to periodic boundary condition
        pos(mdown,d) = pos(mdown,d)-(s*2*pi);
        if isempty(intersect(mup,mdown)) == 0               % if there are values common to both mup and mdown.
            disp('strange');
        end
    end
end

% % % function pos = getRelativePosition(g,v,i)
% % % pos = g.verts(v,1:2);
% % % if(g.bc==1),
% % %     for d=1:2,
% % %         ap = g.verts(v,d);
% % %         [p pid] = min(abs(ap));
% % %         p = ap(pid);
% % %         idx = mod(pid,length(v))+1;
% % %         while (idx~=pid),
% % %             if(abs(pos(idx,d)+2*pi-p)< abs(pos(idx,d)-p)),
% % %                pos(idx,d) =  pos(idx,d)+2*pi;
% % %             end
% % %             if(abs(pos(idx,d)-2*pi-p)< abs(pos(idx,d)-p)),
% % %                pos(idx,d) =  pos(idx,d)-2*pi;
% % %             end
% % %             p= pos(idx,d);
% % %             idx = mod(idx,length(v))+1;
% % %             %[idx pid]
% % %         end
% % %     
% % % %   mup = find(abs(pos(:,d)+2*pi-p)<abs(pos(:,d)-p));
% % % %   mdown = find(abs(pos(:,d)-2*pi-p)<abs(pos(:,d)-p));
% % % %   pos(mup,d) = pos(mup,d)+(2*pi);
% % % %   pos(mdown,d)=pos(mdown,d)-(2*pi);
% % % %   if(length(intersect(mup,mdown))>0)
% % % %       disp('strange');
% % % %   end
% % % 
% % %     end
% % % end


