%FUNCTION FOR PLOTTING DENSITIES
  function M = display_single_fig(nelx, nely, xPhys, loop, P, Pe,...
      Pe_loss, Vavg, je, objhist, volhist, Vbus, Eff, Vbushist)
%   newmap = [1 1 1];
% ncol = size(summer, 1);
% newmap(2:ncol+1,:) = summer;
colormap(gray);

#newmap = summer;
#subplot(2, 2, 1)
#colormap(newmap); 
imagesc(1-xPhys);
#caxis([0 1]); axis equal; axis on;
title(['V_{bus}: ', num2str(sprintf('%0.3f', Vbus)),...
    'V;  Efficiency: ', num2str(sprintf('%4.2f', Eff)), '%']); 

#subplot(2, 2, 4)
#[mx,my]=meshgrid(0:nelx,0:nely);
#imagesc(reshape(Vavg,nely,nelx)); axis equal; axis tight;colorbar; %shading interp;%interpolatie --> shape function;
#subplot(2, 2, 2);
#imagesc(reshape(je,nely,nelx)); axis equal; axis tight;colorbar;drawnow

#subplot(2, 2, 3)
#plot(objhist(1:loop))

#set(gcf,'DoubleBuffer','on'); % double buffer
#set(gcf,'PaperPositionMode','auto');
#M=getframe(gcf);
M = 0;
dummy = 0;
  end