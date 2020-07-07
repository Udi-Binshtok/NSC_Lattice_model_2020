function geo1 = special_MorphologyUpdate(geo,Cells,step,params,n,cells_num,factor)
%SPECIAL_MORPHOLOGYUPDATE will change the morphology of the cells 'cells_num' in the lattice at time step 'step' (current step) 
%   The change in morphology is due to the new cells' states at time step 'step' (that is likely different from time step 'step-1') 
%   The change will be executed by a 2-D vertex model which minimizing the lattice mechanical energy - relaxing it into the new morphology 
    
    geo1 = geo;
    
    %%% scaling areas:
    %%% the area of the lattice and optimal area of each cell type is scaled due to migrated Neurons and new cells from division events
    number_of_cells.initial = length(find(geo(1).g.dead == 0)); % number of cells in the lattice at time step 1 (without deleted cells/migreted Neurons)
    number_of_cells.step = length(find(geo(step).g.dead == 0)); % number of cells in the lattice at time step 'step' (without deleted cells/migreted Neurons)
    geo1(step).g.area_scale = number_of_cells.initial/number_of_cells.step; % scale the lattice total area by the ratio of number of cells: initial/(current step)
    scaled_A0 = (geo1(step).g.area_scale).*(params.A0); % scaling the optimal area for each type of cell state 
    % A0 is a vector of optimal area for each type of cell state [ NSC, aNSC, MC, daughter cells, aNP, Neuron]; 
    
    %%% set optimal areas for each cell according to its state
    NSC = Cells(step).states == 0;              % which cells are NSC at time step 'step' (current step)
    aNSC = Cells(step).states == 1;             % which cells are aNSC at time step 'step' (current step)
    MotherCells = Cells(step).states == 2;      % which cells are Mother cells at time step 'step' (current step)
    DaughterCells = Cells(step).states == 3;    % which cells are daughter cells(Mother cells that has divided) at time step 'step' (current step)
    aNP = Cells(step).states == 4;              % which cells are aNP at time step 'step' (current step)
    div_aNP = Cells(step).states == 4.1;        % which cells are div_aNP (aNP that has divided at time step 'step-1') at time step 'step' (current step)
    Neurons = Cells(step).states == 5;          % which cells are Neurons at time step 'step' (current step)
    Migrated_Neurons = Cells(step).states == -1;% which cells are Migrated_Neurons at time step 'step' (current step)
    
    % factor is a vector consist of scaled_A0 factors for each cell type [NSC, aNSC, MC, DC, delayed_aNP, aNP, div_aNP, Neuron]  
    geo1(step).g.areas(NSC) = factor(1)*scaled_A0(1); 
    geo1(step).g.areas(aNSC) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(MotherCells) = factor(3)*scaled_A0(3);     % mother cells are actually aNSC cells under division process
    geo1(step).g.areas(DaughterCells) = factor(4)*scaled_A0(4);   % daughter cells are actually aNSC cells after division process
    geo1(step).g.areas(aNP) = factor(5)*scaled_A0(5);
    geo1(step).g.areas(div_aNP) = factor(6)*scaled_A0(5);         % div_aNP are aNPs that divide once before differentiating into Neurons
    geo1(step).g.areas(Neurons) = factor(7)*scaled_A0(6);         % because it is out of the lattice
    geo1(step).g.areas(Migrated_Neurons) = 0;
    
    %%% morphlogy changes
    g = geo1(step).g;
    g1 = special_relaxLattice(g,n,params.resolution_sp,cells_num); % update morphology only for the cells 'cells_num'
    
    geo1(step).g = g1;
end

function g = special_relaxLattice(g,n,resolution,cells_num)
% relaxLattice will change the morphology of the lattice by minimizing the mechanical energy of the lattice
%   g is the geometric structure of the lattice
%   n is the number of iterations the function runs
%   resolution is a factor of how fine the change in morphology will be in each iteration
%   cells_num is a vector containing the cells' number on which the function acts 
    
    for i = 1:n
        ve = extractverts(g);
        dE = special_denergy(ve,g,cells_num);
        noE = norm(dE);
        if noE > 1E-10
            dE = resolution*dE/noE;
% %             dE = 0.01*dE/noE;
% %             dE = 0.03*dE/noE;
            ve = ve - dE';
            g = insertverts(ve,g);
        end
    end
end

function dE = special_denergy(v,gin,cells_num)

    dE = zeros(2*length(gin.verts),1);
    g = insertverts(v,gin);
    for i = cells_num
        if g.dead(i) == 0
% %             dEi = dHarea(g,i) + dHlength(g,i) + dHpot(g,i);
            dEi = dHarea(g,i) + dHlength_Roie(g,i);
            dE = dE + dEi;
        end
    end
    dE = dE';
end