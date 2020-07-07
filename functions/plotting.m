function [Movie,MovieName] = plotting(geo,Cells,time)
%PLOTTING will plot and save a movie of the simulation as .avi file in the Matlab current folder 
%   Each cell is at one of these states:
%   NSC (0), aNSC (1), MotherCell (2), DaughterCell (3), aNP (4), div_aNP (4.1), Neuron (5), deleted cell/migrated Neuron (-1)
%   state 4.1 is a div_aNP cell which is an aNP that has divided once and will differenitate on the next step (it will be presented as aNP)
    
    FrameIndex = 1;
    NumOfFrames = 10; % number of frames for each time step;
    MovieName = 'Dray_et_al_2020_SimulationMovie';
    
    edgecolor = [0.5 0.5 0.5];
    
    for t = time
        states = Cells(t).states;
        g = geo(t).g;
        
        cells_in_lattice = find(~g.dead);
        
        F = figure('color','white');
        
        for c = cells_in_lattice'
            verts = (g.bonds(g.cells{c+1}(:),1));
            if(length(verts)>2)
                %%% cell state
                cell_state = states(c);
                switch cell_state
                    case 0
                        facecolor = [ 0 , 1 , 0 ];  % NSC are green
                    case 1
                        facecolor = [ 0.498039 , 0 , 1 ];  % aNSC are violet 
                    case 4
                        facecolor = [ 1 , 0.65 , 0 ];  % aNP are orange
                    case 4.1
                        facecolor = [ 1 , 0.65 , 0 ];  % aNP that has divided once is still aNP, so it is orange
                end

                %%% cell position
                v = getRelativePosition(g,verts,c);

                %%% patching the cell in the lattice (graphically)
                patch(v(:,1),v(:,2),'w','EdgeColor',edgecolor,'FaceColor',facecolor);
                hold on;
                % uncomment to display cell number
%                 text(mean(v(:,1)),mean(v(:,2)),num2str(c),'HorizontalAlignment','center');
                hold on;
            end
        end

        %%% percentage of cells type in the cells lattice
        text_position = F.Children.Position;
        text_position(1) = 1.2*text_position(1); % position x of the bottom-left corner of the text box
        text_position(2) = 0.2*text_position(2) + 0.01*text_position(4); % position y of the bottom-left corner of the text box
        text_position(3) = 1.15*text_position(3); % length of text box on x axis (i.e. figure x axis)
        text_position(4) = 0.08*text_position(4); % length of text box on y axis (i.e. figure y axis)
        
        denominator = length(find(states == 0)) + length(find(states == 1)) + length(find(states == 4)) + length(find(states == 4.1));
        per.aNSC = 100*length(find(states == 1))/denominator;
        per.aNP = 100*(length(find(states == 4)) + length(find(states == 4.1)))/denominator;
        tot_num_cells = length(g.dead(g.dead == 0));
        
        string = {['Time step: ' num2str(t) '  {\color[rgb]{0.498039 0 1}aNSC }' num2str(round(per.aNSC,1)) ' %  {\color[rgb]{1 0.65 0}aNP }' num2str(round(per.aNP,1)) ' %   #Cells: ' num2str(tot_num_cells)]};

        annotation('textbox',text_position,'String',string, 'FontSize',12,'EdgeColor',[1 1 1],'LineWidth',1,'BackgroundColor',[1 1 1])
        
        title(['{\color{green}qNSC }' '{\color[rgb]{0.498039 0 1}aNSC }' '{\color[rgb]{1 0.65 0}aNP }'],'FontSize',14)
            
        axis image; axis off; box off;
        
%         F.WindowState = 'maximized';
        
        %%% creating a movie
        get = getframe(F);
        for  m = FrameIndex:FrameIndex+NumOfFrames-1
            Movie(m) = get; %generates a movie variable
        end
        FrameIndex = FrameIndex + NumOfFrames;
        
        close(F)
        
    end
    v = VideoWriter(MovieName,'Uncompressed AVI');
    open(v)
    writeVideo(v,Movie); % save movie in avi format
    close(v)
end