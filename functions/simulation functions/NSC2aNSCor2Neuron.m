function [geo1,Cells1] = NSC2aNSCor2Neuron(geo,Cells,step,params)
%NSC2ANSC turns NSC at time step 'step-1' to aNSC at time step 'step', with some probability.
%   geo is the geometric strcture (geo(t).g = morphology of the lattice at time step 't'). 
%   Cells is the model structure (Cells(t).states = states of each cell at time step 't'). 
%   step is the current time step.
%   params is the parameter structure.

    Cells1 = Cells;
    geo1 = geo;
    g = geo(step).g;
       
    rate.NSC2aNSC = params.rate_NSC2aNSC;
    rate.rate_NSC2Neuron = params.rate_NSC2Neuron;
    LI = params.LI;                         % Shrot range Notch mediated Lateral Inhibition - strength of inhibition by neighbor aNP
    
    NSC = find(Cells1(step-1).states == 0); % which cells are NSC at time step 'step-1'
    if ~isempty(NSC)
        
        aNSC = zeros(1,length(NSC)); % new aNSC at step 'step' (from the NSCs at time step 'step-1')
        Neurons = zeros(1,length(NSC)); % new Neurons at step 'step' (from the NSCs at time step 'step-1')
        
        %%% regulation by aNPs
        aNP = find(Cells1(step-1).states == 4); % which cells are aNP at time step 'step-1'
        div_aNP1 = find(Cells1(step-1).states == 4.1); % which cells are div_aNP1 (aNP that has divided once) at time step 'step-1'
        aNP = [ aNP div_aNP1 ];
                
        NSC2Neuron = []; % some NSC turn directly into Neuron
        
        for n = 1:length(NSC)
            %%% find neighbors of NSC and check if they are aNP
            Neighbors = g.bonds(g.cells{NSC(n)+1},4);
            aNP_Neighbors = intersect(Neighbors,aNP);

            %%% probability to turn NSC to Neuron (exponential distribution probability density function)
            rate_NSC2Neuron = rate.rate_NSC2Neuron;
            prob = 1 - exp(-rate_NSC2Neuron);
                
            %%% compering probability to random (0,1) num (uniform distribution)
            if prob >= rand
                Neurons(n) = 1; % NSC(n) turns to Neuron
                NSC2Neuron = [NSC2Neuron NSC(n)]; % record NSC->Neuron event
            else
                
                %%% probability to turn NSC to aNSC (exponential distribution probability density function)
                rate_eff = 1/((1/rate.NSC2aNSC) + LI*(length(aNP_Neighbors) ));
                prob = 1 - exp(-rate_eff);

                %%% compering probability to randon num between 0 and 1 (uniform distribution)
                if prob >= rand
                    aNSC(n) = 1; % NSC(n) turns to aNSC
                end
            end
        end

        NewActivatedNSC = NSC(aNSC == 1);
        Cells1(step).states(NewActivatedNSC) = 1;
        NewNeurons = NSC(Neurons == 1);
        Cells1(step).states(NewNeurons) = 5;
        Cells1(step).NSC2Neuron =  NSC2Neuron; % record NSC->Neuron events
        
        mi = params.mi.NSC2aNSCor2Neuron;
%         factor = params.factor;
        factor = [ 1 0.01 1 1 1 1 1 ];
        geo1 = MorphologyUpdate(geo1,Cells1,step,params,mi,factor);
    end
end

