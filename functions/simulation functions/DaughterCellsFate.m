function [geo1,Cells1] = DaughterCellsFate(geo,Cells,step,params)
%DAUGHTERCELLSFATE determines the fate of Daughter cells from time step 'step' 
%   Daughter cells at time step 'step' from Mother cells at time step 'step' can turn into aNSC or aNP 
        
    Cells1 = Cells;
    geo1 = geo;
    g = geo1(step).g;
    
    DC_pairs = Cells1(step).MC_division;
    if ~isempty(DC_pairs)
        
        p_2NSC = params.probability_for_2DaughterCells_to_become_2NSC;
        p_1NSC1aNP = params.probability_for_2DaughterCells_to_become_1NSC1aNP;
        p_2aNP = params.probability_for_2DaughterCells_to_become_2aNP;
        
        for p = 1:size(DC_pairs,1)
            if p_2aNP > rand
                Cells1(step).states(DC_pairs(p,1)) = 4;
                Cells1(step).states(DC_pairs(p,2)) = 4;
            else
                sub_p_tot = p_2NSC + p_1NSC1aNP;
                sub_p_2aNSC = p_2NSC/sub_p_tot;
                if sub_p_2aNSC > rand
                    Cells1(step).states(DC_pairs(p,1)) = 0;
                    Cells1(step).states(DC_pairs(p,2)) = 0;
                else
                    if rand > 0.5
                        Cells1(step).states(DC_pairs(p,1)) = 0;
                        Cells1(step).states(DC_pairs(p,2)) = 4;
                    else
                        Cells1(step).states(DC_pairs(p,1)) = 4;
                        Cells1(step).states(DC_pairs(p,2)) = 0;
                    end
                end
            end
            % special case morphology update - increasing daugher cell size 
            mi_sp = params.mi.DaughterCellsFate_sp;
            factor = params.factor;
            Neighbors_of_cell_DC1 = g.bonds(g.cells{DC_pairs(p,1)+1},4);
            Neighbors_of_cell_DC2 = g.bonds(g.cells{DC_pairs(p,2)+1},4);
            geo1 = special_MorphologyUpdate(geo,Cells1,step,params,mi_sp,union(Neighbors_of_cell_DC1',Neighbors_of_cell_DC2'),factor);

            mi = params.mi.DaughterCellsFate;
            factor = params.factor;           
            geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi,factor);
        end
    end

end

