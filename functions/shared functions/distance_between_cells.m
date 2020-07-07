function [ distance_info ] = distance_between_cells( g_A, g_B , cell_number_A, cell_number_B )
%DISTANCE_BETWEEN_CELLS will calculate the distance between pair of cells within the lattice and also the 'imaged' lattices created from the periodic boundry condition
% 
% cell_number_A is the cell that we measure the distance from it to another cell (to cell_number_B)
% cell_number_A is alwayes in the "original" lattice (from original + imaginary lattices in the periodic boundary condition)
% cell_number_B is once in the "original" lattice and also in the imaginary lattices
% 
% % % % distance_info is a structure containing the information of the distances.
% % % % 
% % % % every cell has a (0;3) tensor: 9 times the (number_of_cell X number_of_cell) matrix.
% % % % 5 is the original lattice.
% % % % 1 is the upper-left imaged lattice.
% % % % 2 is upper-midle. 
% % % % 3 is upper-right.
% % % % 4 is midle-left 
% % % % ... 
% % % % 9 is lower-right
% % % %
% % % % that means that for each cell this function will plot the distance to any
% % % % other cell including cells in imaged lattices
% % % %
% % % % [ cells_info ] = Map_of_cells( g,area );

[ location_A ] = Location_of_one_cell( g_A, cell_number_A );
[ location_B ] = Location_of_one_cell( g_B, cell_number_B ); % NOTE: if g_A is not equal to g_B, the length scale can be different between them so the distance is inaccurate.
                                                             %       See correction for this issue in delta_x and delta_y below    

% go through all lattices (the "original" and the imaginary due to periodic boundary condition)
% (2,2) is the original lattice
% (1,1) is the upper-left imaged lattice, (1,2) is the upper one, and so on
for i = 1:3 
    for j = 1:3
%         delta_x = location_A(2,2,1) - location_B(i,j,1);
%         delta_y = location_A(2,2,2) - location_B(i,j,2);
        delta_x = location_A(2,2,1) - sqrt(g_A.area_scale/g_B.area_scale)*location_B(i,j,1);
        delta_y = location_A(2,2,2) - sqrt(g_A.area_scale/g_B.area_scale)*location_B(i,j,2);
        Dist = sqrt(delta_x^2 + delta_y^2);
        
        distance_info(i,j,1) = delta_x; % distance on x axis
        distance_info(i,j,2) = delta_y; % distance on y axis
        distance_info(i,j,3) = Dist; % distance between the cells (center to center)
    end
end
    
% % % k = size(cells_info,2); % k is the # of cells
% % % Dist = NaN(k,k,9);
% % % 
% % % % I inserted these two arrays for killing filopodia in the function fixed filopodia
% % % v_x = NaN(k,k,9);
% % % v_y = NaN(k,k,9);
% % % 
% % % for i = 1:k    
% % %     for j = 1:k
% % %         % go through all lattices
% % %         Dist(i,j,1) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,1,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,1,2))^2);
% % %         v_x(i,j,1) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,1,1);
% % %         v_y(i,j,1) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,1,2);
% % %         
% % %         Dist(i,j,2) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,2,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,2,2))^2);
% % %         v_x(i,j,2) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,2,1);
% % %         v_y(i,j,2) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,2,2);
% % %         
% % %         Dist(i,j,3) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,3,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,3,2))^2);
% % %         v_x(i,j,3) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(1,3,1);
% % %         v_y(i,j,3) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(1,3,2);
% % %         
% % %         Dist(i,j,4) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,1,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,1,2))^2);
% % %         v_x(i,j,4) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,1,1);
% % %         v_y(i,j,4) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,1,2);
% % %         
% % %         Dist(i,j,5) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,2,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,2,2))^2);
% % %         v_x(i,j,5) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,2,1);
% % %         v_y(i,j,5) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,2,2);
% % %         
% % %         Dist(i,j,6) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,3,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,3,2))^2);
% % %         v_x(i,j,6) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(2,3,1);
% % %         v_y(i,j,6) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(2,3,2);
% % %         
% % %         Dist(i,j,7) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,1,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,1,2))^2);
% % %         v_x(i,j,7) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,1,1);
% % %         v_y(i,j,7) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,1,2);
% % %         
% % %         Dist(i,j,8) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,2,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,2,2))^2);
% % %         v_x(i,j,8) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,2,1);
% % %         v_y(i,j,8) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,2,2);
% % %         
% % %         Dist(i,j,9) = sqrt((cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,3,1))^2 + (cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,3,2))^2);
% % %         v_x(i,j,9) = cells_info(1,i).position(2,2,1) - cells_info(1,j).position(3,3,1);
% % %         v_y(i,j,9) = cells_info(1,i).position(2,2,2) - cells_info(1,j).position(3,3,2);
% % %         
% % %     end
% % %     % distance_info{1,i} = Dist(i,:,:); % that is not needed
% % % end
% % % 
% % % distance_info = Dist;
% % % distance_values_x = v_x;
% % % distance_values_y = v_y;

end

