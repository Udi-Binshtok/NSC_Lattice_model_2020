function [params] = defaultparams(g)
    
    %%% simulation number of time steps
    params.Time_max = 500; % max time step
    
    %%% Initial conditions parameters
    params.per_aNSC = 3.08; % observed percentage of aNSC in the Medial pallium (DM)
    params.per_aNP = 12; %20; %16; % observed percentage of aNP in the Medial pallium (DM)
        
    %%% short range regulation
    params.LI = 1000; % Shrot range Notch mediated Lateral Inhibition - strength of inhibition by neighbor aNP

    %%% rates
    params.rate_NSC2aNSC =          0.014;  % qNSC activation rate
    params.rate_NSC2Neuron =        0.0014; % qNSC direct differentiation rate
    params.rate_aNSC2MotherCell =   0.22;   % aNSC division rate
    params.rate_aNP2Neuron =        0.062;  % aNP differentiation rate 
    params.rate_aNP_division =      0.032;  % aNP division rate
    params.rate_divaNP2Neuron =     0.062;  % div_aNP differentiation rate 
    
    %%% probabilities
    params.probability_for_2DaughterCells_to_become_2NSC = 0.30;     % Symmetric division 2 qNSC
    params.probability_for_2DaughterCells_to_become_1NSC1aNP = 0.58; % Asymmetric division 1 qNSC 1 aNP
    params.probability_for_2DaughterCells_to_become_2aNP = 0.12;     % Symmetric division 2 aNP
    
    
    %% Morphology update
    
    %%% set g.paras corresponding to the coefficents of area, perimeter and bonds' length in the mechanical energy minimization function  
    params.g_paras = [ 40 ; 2.5 ; 5 ]; % [ area ; perimeter ; bonds' length ] 
    
    %%% set the optimal area for each cell type (this area is scaled in each time step according to the number of cells in the lattice at that time step relative to the initial number of cells)                               
    tot_area = (2*pi)*(2*pi); % the length and height of the lattice each equal 2pi 
    optimal_area.NSC = tot_area/length(find(g.dead == 0));
    optimal_area.aNSC =             0.35*optimal_area.NSC;
    optimal_area.MC =               0.35*optimal_area.NSC; 
    optimal_area.DaughterCells =    0.175*optimal_area.NSC; 
    optimal_area.aNP =              0.1*optimal_area.NSC; 
    optimal_area.Neuron =           0.05*optimal_area.NSC; 
    params.A0 = [ optimal_area.NSC , optimal_area.aNSC , optimal_area.MC , optimal_area.DaughterCells,  optimal_area.aNP , optimal_area.Neuron ];
    % see also factor in MorphologyUpdate and special_MorphologyUpdate (increase or decrease in optimal area, by a factor) 
    params.factor = [ 1 1 1 1 1 1 1 ];  % factor is a vector consist of optimal area multiplying factors for each cell type [NSC, aNSC, MC, DC, aNP, div_aNP, Neuron]
                   
    %%% Morphology update DeltaE resolution (see relaxLattice function and also special_relaxLattice within special_MorphologyUpdate)
    params.resolution = 0.001;
    params.resolution_sp = 0.001; 
    
    %%% number of morphology update iterations (mi) for each function in the simulation function
    % sp menas special case where morphology changes only aroound selected cells
    global_itreations = 0;
    mi.get_IC =                         200;
    mi.NSC2aNSCor2Neuron =              200 + global_itreations;
    mi.aNSC2MotherCell =                10 + global_itreations;
    mi.MotherCells2DaughterCells_sp =   50 + global_itreations; 
    mi.MotherCells2DaughterCells =      5 + global_itreations;
    mi.DaughterCellsFate_sp =           200 + global_itreations; 
    mi.DaughterCellsFate =              10 + global_itreations;
    mi.aNPsFate1_divide_sp =            200 + global_itreations;
    mi.div_aNPs1_sp =                   50 + global_itreations;
    mi.aNPsFate1_sp =                   200 + global_itreations;
    mi.aNPsFate1 =                      10 + global_itreations;
    mi.NeuronMigation =                 90 + global_itreations; 
    mi.simulation =                     50 + global_itreations;
    
    params.mi = mi;
    
end

