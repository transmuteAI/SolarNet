function freelems = get_freelems(design, edofMat, actualdofs, nely, nelx)
    freelems = [];
    if design == -1 || design == 0
        freelems = setdiff(1:nely*nelx, freelems);
    else
%         for i = 1:1:nelx*nely
%             if(sum(ismember(edofMat(i, :), actualdofs)) == 4)
%                 freelems = [freelems i];
%                 i
%             end
%         end
%         save freelems_hexcell1.mat freelems;

        %load freelems.mat;
        %load square_pum.mat;
        %load freelems_maple.mat;
        %load freelems_hexagon.mat;
        %load freelems_hexagon_6node.mat;
        %load freelems_bike_head.mat;k
        %load freelems_bike_head2.mat;
        % load freelems_bike_side.mat;
        %load freelems_bike_side2.mat
        %load freelems_circle.dat
        
        %load freelems_christmastree.mat; % This one was for the 2 times
        %load freelems_christmastree1.mat; % 1 time without symmetry
        %load freelems_christmastree2.mat;
        load freelems_hexcell1.mat;
        %higher resolution
        
        
    end
end