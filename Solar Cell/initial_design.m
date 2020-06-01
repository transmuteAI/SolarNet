%function for the initial design
function x = initial_design(nelx, nely, elec, volfrac, evec)

global evec2d;

if elec == 1
    x = zeros(nely,nelx);
   x(:, :) = volfrac;
   selected_points = round(linspace(1, 299, 150));
%x(evec,495:515) = 1;
x(evec, 1:5) = 1;
%for the leaf
% x(30:859, 210:225) = 1;
% 
% center_right = 219;
% center_left = 219;
% yline = [180 300 420 540 660 780];
% times = 30;
% deltax = 5;
% deltay = 5;
% for i = 1:1:times
%     if i == 13 
%         yline = yline([2 3 4 5 6]);
%     elseif i == 20
%         yline = yline([2 3 4]);
%     elseif i == 25
%         yline = yline([2 3]);
%     elseif i == 28
%         yline = yline([1]);
%     end
%     for i = 1:1:deltay
%         x(yline+(i-1), center_right:center_right + deltax) = 1;
%         x(yline+(i-1), center_left - deltax:center_left) = 1;
%     end
%     yline = yline - (deltay - 1);
%     center_right = center_right + deltax;
%     center_left = center_left - deltax;
% end

%(805:810, 170:270) = 1;
else
x = repmat(volfrac,nely,nelx);
end


end