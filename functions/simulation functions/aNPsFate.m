function [geo1,Cells1,StopSimulation] = aNPsFate(geo,Cells,step,params)
%ANPSFATE determines the fate of aNP and div_aNP from time step 'step-1'
%   aNP can divide once to generate 2 aNPs that will differentiate into 2 neurones or directly differentiate into neuron
   
    StopSimulation = 0;
        
    geo1 = geo;
    g = geo1(step).g;
    Cells1 = Cells;
    
    %%% aNP that did not divide or differentiated yet
    aNP = find(Cells1(step-1).states == 4); % which cells are aNP at time step 'step-1'
    if ~isempty(aNP)
       
        for n = 1:length(aNP)
            %%% probability for aNP to divide (exponential distribution probability density function)
            rate_aNPdivision = params.rate_aNP_division;
            prob = 1-exp(-rate_aNPdivision);
            
            if prob >= rand
                % aNP division
                mi_sp1 = params.mi.aNPsFate1_divide_sp;
                factor = params.factor;
                factor(5) = 20; % %NOTE: Increasing aNP area before division. factor is a vector consist of prefarable area factors for each cell type [NSC, aNSC, MC, DC, aNP, div_aNP, Neuron]
                Neighbors_of_cell_n = g.bonds(g.cells{aNP(n)+1},4);
                geo1 = special_MorphologyUpdate(geo1,Cells1,step,params,mi_sp1,[aNP(n) Neighbors_of_cell_n'],factor); % special_MorphologyUpdate changes the morphology only to the cells indicated at the 6th entrance of the function (see the function description)
                
                [g, new_cell_num, weird] = cell_division( g, aNP(n) ); % divide aNP
                if weird
                    disp('weird - some cells have only 2 bonds')
                    StopSimulation = 1;
                    return
%                     keyboard
                end
                Cells1(step).states(aNP(n)) = 4.1;              % state of all aNP from time step 'step-1' that divide is '4.1' at time step 'step'
                Cells1(step).states(new_cell_num) = 4.1;        % state of all new cells (aNP) is '4.1'
                g.dead(new_cell_num) = 0;                       % new cells are not dead
                g.xboundary(new_cell_num,1:2) = 0;              % 0 means periodic boundary condition
                g.yboundary(new_cell_num,1:2) = 0;              % 0 means periodic boundary condition
                geo1(step).g = g;
                
                mi_sp2 = params.mi.div_aNPs1_sp;
                factor = params.factor;
                factor(6) = 20; % NOTE: Increasing divided_aNP area after division. factor is a vector consist of prefarable area factors for each cell type [NSC, aNSC, MC, DC, aNP, div_aNP, Neuron]
                Neighbors_of_cell_n = g.bonds(g.cells{aNP(n)+1},4);
                geo1 = special_MorphologyUpdate(geo1,Cells1,step,params,mi_sp2,[aNP(n) new_cell_num Neighbors_of_cell_n'],factor); % special_MorphologyUpdate changes the morphology only to the cells indicated at the 6th entrance of the function (see the function description)
                
            else
                % if the aNP did not divide than it can differentiates into a Neuron
                %%% probability to turn aNP to Neuron (exponential distribution probability density function)
                rate_aNP2Neuron = params.rate_aNP2Neuron;
                prob = 1 - exp(-rate_aNP2Neuron);
                
                %%% compering probability to randon num between 0 and 1 (uniform distribution)
                if prob >= rand
                    Cells1(step).states(aNP(n)) = 5;
                    
                    factor = params.factor;
                    mi_sp_aNP2N = params.mi.aNPsFate1_sp;
                    Neighbors_of_cell_n = g.bonds(g.cells{aNP(n)+1},4);
                    geo1 = special_MorphologyUpdate(geo1,Cells1,step,params,mi_sp_aNP2N,[aNP(n) Neighbors_of_cell_n'],factor); % special_MorphologyUpdate changes the morphology only to the cells indicated at the 6th entrance of the function (see the function description)
                end
            end
            
            mi1 = params.mi.aNPsFate1;
            factor = params.factor;
            geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi1,factor);
        end
    end
    
%%
    %%% aNP that already divided once can differentiate into Neuron
    div_aNP = find(Cells1(step-1).states == 4.1); % which aNP cells has already divided and at state 4.1 at time step 'step-1'
    if ~isempty(div_aNP) 
        for dn = 1:length(div_aNP) 
            %%% probability to turn aNP to Neuron (exponential distribution probability density function)
            rate_divaNP2Neuron = params.rate_divaNP2Neuron;
            prob = 1-exp(-rate_divaNP2Neuron);

            %%% compering probability to randon num between 0 and 1 (uniform distribution)
            if prob >= rand
                Cells1(step).states(div_aNP(dn)) = 5;

                mi1 = params.mi.aNPsFate1;
                factor = params.factor;
                geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi1,factor);
            end
            
            mi1 = params.mi.aNPsFate1;
            factor = params.factor;
            geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi1,factor);
        end
    end
 
end

