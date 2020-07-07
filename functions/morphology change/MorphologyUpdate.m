function geo1 = MorphologyUpdate(geo,Cells,step,params,n,factor)
%MORPHOLOGYUPDATE will change the morphology of the lattice at time step 'step' (current step) 
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
    
    % factor is a vector consist of scaled_A0 factors for each cell type [NSC, aNSC, MC, DC, aNP, div_aNP, Neuron]
    
    % scale the NSC factor by the area of the cell (smaller cells will get bigger A0 and larger cells get smaller A0)
    NSC_cell_num = find(NSC == 1);
    NSC_area = NaN(1,length(NSC_cell_num));  
    for i = 1:length(NSC_cell_num)
        NSC_area(i) = cellarea(geo1(step).g,NSC_cell_num(i));
    end
    NSC_av_area = mean(NSC_area);
    NSC_expo = 3;
    for i = 1:length(NSC_cell_num)
        geo1(step).g.areas(NSC_cell_num(i)) = ((NSC_av_area/NSC_area(i))^NSC_expo)*factor(1)*scaled_A0(1);
    end 
    
    % scale the aNSC factor by the area of the cell (smaller cells will get bigger A0 and larger cells get smaller A0)
    aNSC_cell_num = find(aNSC == 1);
    aNSC_expo = 4;
    for i = 1:length(aNSC_cell_num)
        aNSC_area = cellarea(geo1(step).g,aNSC_cell_num(i));
        geo1(step).g.areas(aNSC_cell_num(i)) = ((scaled_A0(2)/aNSC_area)^aNSC_expo)*factor(2)*scaled_A0(2);
    end
%     aNSC_cell_num = find(aNSC == 1);
%     aNSC_area = NaN(1,length(aNSC_cell_num));  
%     for i = 1:length(aNSC_cell_num)
%         aNSC_area(i) = cellarea(geo1(step).g,aNSC_cell_num(i));
%     end
%     aNSC_av_area = mean(aNSC_area);
%     aNSC_expo = 4;
%     for i = 1:length(aNSC_cell_num)
%         geo1(step).g.areas(aNSC_cell_num(i)) = ((aNSC_av_area/aNSC_area(i))^aNSC_expo)*factor(2)*scaled_A0(2);
%     end
    
    aNP_cell_num = find(aNP == 1);
    aNP_expo = 4;
    for i = 1:length(aNP_cell_num)
        aNP_area = cellarea(geo1(step).g,aNP_cell_num(i));
        geo1(step).g.areas(aNP_cell_num(i)) = ((scaled_A0(5)/aNP_area)^aNP_expo)*factor(5)*scaled_A0(5);
    end
    
    div_aNP_cell_num = find(div_aNP == 1);
    div_aNP_expo = 4;
    for i = 1:length(div_aNP_cell_num)
        div_aNP_area = cellarea(geo1(step).g,div_aNP_cell_num(i));
        geo1(step).g.areas(div_aNP_cell_num(i)) = ((scaled_A0(5)/div_aNP_area)^div_aNP_expo)*factor(6)*scaled_A0(5);
    end

%     geo1(step).g.areas(NSC) = factor(1)*scaled_A0(1);
%     geo1(step).g.areas(aNSC) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(MotherCells) = factor(3)*scaled_A0(3);     % mother cells are actually aNSC cells under division process
    geo1(step).g.areas(DaughterCells) = factor(4)*scaled_A0(4);   % daughter cells are actually aNSC cells after division process
%     geo1(step).g.areas(aNP) = factor(5)*scaled_A0(5);
%     geo1(step).g.areas(div_aNP) = factor(6)*scaled_A0(5);         % div_aNP are aNPs that divide once before differentiating into Neurons
    geo1(step).g.areas(Neurons) = factor(7)*scaled_A0(6);         % is the area of Neuron actully 0? because it is out of the lattice
    geo1(step).g.areas(Migrated_Neurons) = 0;
    
    %%% morphlogy changes
    g = geo1(step).g;
    g1 = relaxLattice(g,n,params.resolution); % update morphology
    
    geo1(step).g = g1;

end

