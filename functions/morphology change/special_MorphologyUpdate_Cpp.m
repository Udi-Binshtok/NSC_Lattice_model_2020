function geo1 = special_MorphologyUpdate_Cpp(geo,Cells,step,params,n,cells_num,factor)
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
    div_aNSC11 = Cells(step).states == 1.11;    % which cells are div_aNSC at time step 'step'
    div_aNSC12 = Cells(step).states == 1.12;    % which cells are div_aNSC at time step 'step'
    div_aNSC2 = Cells(step).states == 1.2;      % which cells are div_aNSC2 at time step 'step'
    MotherCells = Cells(step).states == 2;      % which cells are Mother cells at time step 'step' (current step)
    MotherCells2 = Cells(step).states == 2.1;   % which cells are Mother cells2 at time step 'step' (current step)
    DaughterCells = Cells(step).states == 3;    % which cells are daughter cells(Mother cells that has divided) at time step 'step' (current step)
    DaughterCells2 = Cells(step).states == 3.2; % which cells are daughter cells2(Mother cells2 that has divided) at time step 'step' (current step)
    delayed_aNP = Cells(step).states == 3.1;    % which cells are delayed_aNP(aNSC that will become aNP) at time step 'step' (current step)
    aNP = Cells(step).states == 4;              % which cells are aNP at time step 'step' (current step)
    div_aNP = Cells(step).states == 4.1;        % which cells are div_aNP (aNP that has divided at time step 'step-1') at time step 'step' (current step)
    Neurons = Cells(step).states == 5;          % which cells are Neurons at time step 'step' (current step)
    Migrated_Neurons = Cells(step).states == -1;% which cells are Migrated_Neurons at time step 'step' (current step)
    
    % factor is a vector consist of scaled_A0 factors for each cell type [NSC, aNSC, MC, DC, delayed_aNP, aNP, div_aNP, Neuron]  
    geo1(step).g.areas(NSC) = factor(1)*scaled_A0(1); 
    geo1(step).g.areas(aNSC) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(div_aNSC11) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(div_aNSC12) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(div_aNSC2) = factor(2)*scaled_A0(2);
    geo1(step).g.areas(MotherCells) = factor(3)*scaled_A0(3);     % mother cells are actually aNSC cells under division process
    geo1(step).g.areas(MotherCells2) = factor(3)*scaled_A0(3);
    geo1(step).g.areas(DaughterCells) = factor(4)*scaled_A0(4);   % daughter cells are actually aNSC cells after division process
    geo1(step).g.areas(DaughterCells2) = factor(4)*scaled_A0(4);
    geo1(step).g.areas(delayed_aNP) = factor(5)*scaled_A0(2);     % delayed_aNP is an aNSC that will become aNP (it is after cell division event)
    geo1(step).g.areas(aNP) = factor(6)*scaled_A0(5);
    geo1(step).g.areas(div_aNP) = factor(7)*scaled_A0(5);         % div_aNP are aNPs that divide once before differentiating into Neurons
    geo1(step).g.areas(Neurons) = factor(8)*scaled_A0(6);         % is the area of Neuron actully 0? because it is out of the lattice
    geo1(step).g.areas(Migrated_Neurons) = 0;
    
    %%% morphlogy changes
    g = geo1(step).g;
    original_paras = g.paras;
    Cpp_paras = params.g_paras;
    g.paras = Cpp_paras;
    %%% send "special_relaxLattice" to run in C++
    %%% g1 = special_relaxLattice(g,n,params.resolution_sp,cells_num); % update morphology only for the cells 'cells_num'
    %% Micha modifications
    addpath('C:\Udi\LaureModel_Cpp\CppFiles'); %% update according to your directory structure
%     addpath('D:\Students\Udi\LaureModel_Cpp\CppFiles'); %% at David lab
    cc=cellfun(@(x) length(x),g.cells);
    nbonds = cc(2:end);
    maxbonds = max(nbonds);
    ncells = length(g.cells);
    cells = -ones(ncells-1,maxbonds);
    for i=1:ncells-1,
        cells(i,1:length(g.cells{i+1})) = g.cells{i+1}-1;%% -1 for 0-indexing (for C++)
    end
%     frozen_cells = find(g.dead == 0);
%     frozen_cells = frozen_cells(~ismember(frozen_cells,cells_num)) - 1;
    frozen_cells = ones(length(Cells(step).states),1); % if a cell is 1 then it is frozen - its vertices stay fixed  
    frozen_cells(cells_num,1) = 0; % if a cell is 0 then it is not frozen and can be changed during energy minimization 
    bonds = g.bonds(:,1:4)-1; %% -1 for 0-indexing (for C++)
%     dur = n/1000;  
    dur = n*params.resolution_Cpp;
%     [v1, b1, c1]= runMechanicalForcesUdi(cells',bonds',nbonds',[dur g.paras'],g.verts(:,1:2)',g.areas',g.bc);
    [v1, b1, c1]= runMechanicalForcesUdi(cells',bonds',nbonds',[dur g.paras'],g.verts(:,1:2)',g.areas',g.bc,frozen_cells);  
    %%%% rebuild g structure from arrays returned by the C++ routine 
        g=insertverts(v1,g);
        c1=c1';
        g.bonds(:,1:4)=b1'+1;
        for i=1:ncells-1,
            g.cells{i+1} = c1(i,c1(i,:)>=0)+1; %% -1 for 0-indexing
        end
    %%%% 	
    g.paras = original_paras;
    geo1(step).g = g;
%     LatticePresentation(g,Cells(step))
end

% % function g = special_relaxLattice(g,n,resolution,cells_num)
% % % relaxLattice will change the morphology of the lattice by minimizing the mechanical energy of the lattice
% % %   g is the geometric structure of the lattice
% % %   n is the number of iterations the function runs
% % %   resolution is a factor of how fine the change in morphology will be in each iteration
% % %   cells_num is a vector containing the cells' number on which the function acts 
% %     
% %     for i = 1:n
% %         ve = extractverts(g);
% %         dE = special_denergy(ve,g,cells_num);
% %         noE = norm(dE);
% %         if noE > 1E-10
% %             dE = resolution*dE/noE;
% % % %             dE = 0.01*dE/noE;
% % % %             dE = 0.03*dE/noE;
% %             ve = ve - dE';
% %             g = insertverts(ve,g);
% %         end
% %     end
% % end
% % 
% % function dE = special_denergy(v,gin,cells_num)
% % 
% %     dE = zeros(2*length(gin.verts),1);
% %     g = insertverts(v,gin);
% %     for i = cells_num
% %         if g.dead(i) == 0
% % % %             dEi = dHarea(g,i) + dHlength(g,i) + dHpot(g,i);
% %             dEi = dHarea(g,i) + dHlength_Roie(g,i);
% %             dE = dE + dEi;
% %         end
% %     end
% %     dE = dE';
% % end