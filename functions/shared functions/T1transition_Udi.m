function g1 = T1transition_Udi(g,b,pert)
    
    %%% see my changes to micha's code at the part: "update vertex position" 

    % b is the bond of cell c (see kill_cells function) that goes T1 transition
    % binv (or b') is the same bond as b but only from the direction of the cell that sahring bond b with cell c, lets call it cell c'
    if (nargin < 3)
        pert = 0; 
    end
    g1 = g;
    binv = find(g.bonds(:,1) == g.bonds(b,2) & g.bonds(:,2) == g.bonds(b,1)); % this is the bond number of binv (or b'), which is the same bond b, but from the side of the neighbor cell c' that shares this bond with cell c
    if length(binv) ~= 1 % it should be only one responding "opposite" bond to bond b
        %  disp('returning');
        return;
    end
    edges = [b binv]; % same bond, once from cell c side (b) and once from its neighbor side c' (binv)
    for v = 1:2
        bi = find(g.bonds(:,v) == g.bonds(b,v)); % find all bonds that shares the same first (v=1) or second (v=2) vertex as b
        bi = setdiff(bi,b); % exclude bond b from bi
        if length(bi) ~= 2 % if the number of bonds, without bond b, that start at vertex g.bonds(b,v) is not 2
            % disp('not returning');
            %  return;
        end
    end
    c1 = g.bonds(b,3); % cell number of cell c
    c3 = g.bonds(b,4); % cell number of cell c', which is cell's c neighbor sharing bond b
%   if(c1*c3==0),
%       disp('not returning');
%       %  return;
%   end
    f1 = find(g.cells{c1+1} == edges(1)); % which of the bonds of cell c is bond b
    if f1 == 1 
        b1 = g.cells{c1+1}(length(g.cells{c1+1})); % choosing the cell's c adjacent bond to bond b
    else
        b1 = g.cells{c1+1}(f1-1); % choosing the cell's c adjacent bond to bond b
    end

    f2 = find(g.cells{c3+1} == edges(2)); % which of the bonds of cell c' is bond binv
    if f2 == 1
        b2 = g.cells{c3+1}(length(g.cells{c3+1})); % choosing the cell's c' adjacent bond to bond binv
    else
        b2 = g.cells{c3+1}(f2-1); % choosing the cell's c' adjacent bond to bond binv
    end
    c2 = g.bonds(b1,4); % cell number of cell's c neighbor sharing bond b1 (adjacent to bond b)
    c4 = g.bonds(b2,4); % cell number of cell's c' neighbor sharing bond b2 (adjacent to bond binv)
    if c2+c4 == 0 || c1 == c4
 %      disp('not returning');
% % %         keyboard
        g1 = ReduceBond(g1,b);
        return;
    end
    if c1*c3 == 0 && c2*c4 == 0
        keyboard
        return;
    end
    
%   if(g.level(c4,1)*g.level(c2,1) > 1),
%   %   disp('returning');
%       return;
%   end   
% %     disp(['flipping ' num2str(b)]);  % uncomment to display 

    %% removing edge from cell c1 and c3
    %g1.cells{c1+1}=setdiff(g.cells{c1+1},edges(1));
    %g1.cells{c3+1}=setdiff(g.cells{c3+1},edges(2));
    % do this instead to avoid sorting and potentially spare some time
    g1.cells{c1+1} = g.cells{c1+1}(~ismember(g.cells{c1+1},edges(1))); % remove bond b from bonds of cell c
    g1.cells{c3+1} = g.cells{c3+1}(~ismember(g.cells{c3+1},edges(2))); % remove bond binv from bonds of cell c'
 
    %% add edge to cell c2 and c4 between edges bidx
    %g.bonds(b1,:)
    b1v = find(g.bonds(:,1) == g.bonds(b1,2) & g.bonds(:,2)== g.bonds(b1,1)); % the number of the "opposite" bond of bond b1 (b1inv or b1')
    bic = find(g.cells{c2+1} == b1v); % find the position of bond b1v (b1') in the bonds of c2
    g1.cells{c2+1} = [g.cells{c2+1}(1:bic-1) edges(1) g.cells{c2+1}(bic:length(g.cells{c2+1}))]; % add bond b to the bonds of c2, at the position of bond b1v (b1')
    b2v = find(g.bonds(:,1) == g.bonds(b2,2) & g.bonds(:,2)== g.bonds(b2,1)); % the number of the "opposite" bond of bond b2 (b2inv or b2')
    bic = find(g.cells{c4+1} == b2v); % find the position of bond b2v (b2') in the bonds of c4
    g1.cells{c4+1} = [g.cells{c4+1}(1:bic-1) edges(2) g.cells{c4+1}(bic:length(g.cells{c4+1}))]; % add bond binv (binv or b') to the bonds of c4, at the position of bond b2v (b2')

    %% move edge extremities
    g1.bonds(b1,2) = g.bonds(edges(1),2);
    g1.bonds(b2,2) = g.bonds(edges(1),1);
    g1.bonds(b1v,1) = g.bonds(edges(2),1);
    g1.bonds(b2v,1) = g.bonds(edges(2),2);

    %% update edge neighbors
    g1.bonds(edges(1),3) = c2;
    g1.bonds(edges(1),4) = c4;
    g1.bonds(edges(2),3) = c4;
    g1.bonds(edges(2),4) = c2;

    %% update vertex position    
% % %     % vert = getRelativePosition(g,g.bonds(b,1:2),c1);
% % %     %mid = 0.5*(vert(1,:)+vert(2,:));
% % %     % g1.verts(g.bonds(b,1),1:2)= mid;
% % %     % g1.verts(g.bonds(b,2),1:2)= mid;
% % %     [blen, mid] = getBoundaryLength(g,b);
% % %     pvec1 = [rand() rand()];
% % %     pvec2 = [rand() rand()];
% % %     pvec1 = pvec1/norm(pvec1);
% % %     pvec2 = pvec2/norm(pvec2);
% % %     g1.verts(g.bonds(b,1),1:2) = mid+pert*blen*pvec1;
% % %     g1.verts(g.bonds(b,2),1:2) = mid+pert*blen*pvec2;
% % % 
% % %     %[c1 c2 c3 c4]
% % %     %g1.cells{c2+1} = sort_edges(g1,c2+1);
% % %     %g1.cells{c4+1} = sort_edges(g1,c4+1);
% % %     %g1.cells{c1+1} = sort_edges(g1,c1+1);
% % %     %g1.cells{c3+1} = sort_edges(g1,c3+1);
% % %     %[c1 c2 c3 c4]

    %%% this is my way of doing this part
    % the length (blen) and center position (mid) of bond b (and also of bond b' or binv)
    [blen, mid] = getBoundaryLength(g,b); 
    
    % get the geometrical center of new cell c1 (that is, cell c1 without bond b)
    c1_new_bonds = g1.cells{c1+1};                                               
    c1_new_verts_indx = g1.bonds(c1_new_bonds,1);                                      
    c1_new_verts_coordinates = getRelativePosition(g1,c1_new_verts_indx,c1);
    c1_new_CM = mean(c1_new_verts_coordinates); 
    vector_mid_to_c1 = c1_new_CM - mid; % this is the vector from mid to c1 center
    % take the position of the second vertex of bond b towards the center of cell c1 
    % pert will be the amount of movment along the mid-c1_center axis
    g1.verts(g.bonds(b,2),1:2) = mid + pert.*vector_mid_to_c1;
    
    % get the geometrical center of new cell c3 (that is, cell c3 without bond b' or binv)
    c3_new_bonds = g1.cells{c3+1};                                               
    c3_new_verts_indx = g1.bonds(c3_new_bonds,1);                                      
    c3_new_verts_coordinates = getRelativePosition(g1,c3_new_verts_indx,c3);
    c3_new_CM = mean(c3_new_verts_coordinates); 
    vector_mid_to_c3 = c3_new_CM - mid; % this is the vector from mid to c3 center
    % take the position of the first vertex of bond b towards the center of cell c3
    % pert will be the amount of movment along the mid-c3_center axis
    g1.verts(g.bonds(b,1),1:2) = mid + pert.*vector_mid_to_c3;
    

    %% update boundary information
    if(g.bc)
        for i = 1:2
            if g.xboundary(c2,i) && g.xboundary(c1,3-i) && g.xboundary(c3,3-i)
                g1.xboundary(c4,3-i) = 1; 
            end
            if g.yboundary(c2,i) && g.yboundary(c1,3-i) && g.yboundary(c3,3-i)
                g1.yboundary(c4,3-i) = 1; 
            end

            if g.xboundary(c4,i) && g.xboundary(c1,3-i) && g.xboundary(c3,3-i)
                g1.xboundary(c2,3-i) = 1; 
            end
            if g.yboundary(c4,i) && g.yboundary(c1,3-i) && g.yboundary(c3,3-i)
                g1.yboundary(c2,3-i) = 1; 
            end
            if g.xboundary(c1,i) && g.xboundary(c2,i)&& g.xboundary(c4,i) && g.xboundary(c3,3-i)
                g1.xboundary(c1,i) = 0; 
            end
            if g.yboundary(c1,i) && g.yboundary(c2,i)&& g.yboundary(c4,i) && g.yboundary(c3,3-i)
                g1.yboundary(c1,i) = 0; 
            end
            if g.xboundary(c3,i) && g.xboundary(c2,i)&& g.xboundary(c4,i) && g.xboundary(c1,3-i)
                g1.xboundary(c3,i) = 0; 
            end
            if g.yboundary(c3,i) && g.yboundary(c2,i)&& g.yboundary(c4,i) && g.yboundary(c1,3-i)
                g1.yboundary(c3,i) = 0; 
            end
        end
    end
end
