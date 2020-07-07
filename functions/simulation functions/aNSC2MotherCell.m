function [geo1,Cells1] = aNSC2MotherCell(geo,Cells,step,params)
%ANSC2MOTHERCELLOR2NSC turns aNSC at time step 'step-1' to Mother cell or NSC at time step 'step', with some probability.
%   geo is the geometric strcture (geo(t).g = morphology of the lattice at time step 't').
%   Cells is the model structure (Cells(t).states = states of each cell at time step 't'). 
%   step is the current time step.
%   params is the parameter structure.
    
    Cells1 = Cells;
    geo1 = geo;
    
    aNSC = find(Cells1(step-1).states == 1);    % which cells are aNSC at time step 'step-1'
    if ~isempty(aNSC)
    
        MotherCells = zeros(1,length(aNSC));    % new MotherCells at step 'step' (from the aNSCs at time step 'step-1')
               
        for a = 1:length(aNSC)
            %%% probability to turn aNSC to Mother Cell (exponential distribution probability density function)
            rate_aNSC2MotherCell = params.rate_aNSC2MotherCell;
            prob = 1 - exp(-rate_aNSC2MotherCell);
            
            %%% comparing probability to random (0,1) num (uniform distribution)
            if prob >= rand
                MotherCells(a) = 1; % aNSC(i) turns to Mother Cell
            end   
        end

        NewMotherCells = aNSC(MotherCells == 1);
        Cells1(step).states(NewMotherCells) = 2;
        
        mi = params.mi.aNSC2MotherCell;
        factor = params.factor;
        geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi,factor);
    end
    
end

