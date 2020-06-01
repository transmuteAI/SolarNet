%FUNCTION FOR PLOTTING DENSITIES
  function M = plot_densities(nelx, nely, xPhys, loop, P, Pe,...
      Pe_loss, Vavg, je, objhist, volhist, Vbus, Eff, Vbushist)
  newmap = [1 1 1];
% ncol = size(copper, 1);
% newmap(2:ncol+1,:) = copper;
% %newmap = summer;
newmap = copper;
ss=get(0,'screensize');
figure(1); %set(1,'position',[6 2*ss(4)/4+20 ss(3)/3 ss(4)/3]);

colormap('gray'); imagesc(1-xPhys);
rectangle('POsition', [0.365 0.82 nelx+0.2 nely-0.32], 'edgecolor', [0.7 0.6 0.5], 'LineWidth', 2);
caxis([0 1]); axis equal; axis on;
%set(gca,'position',[0 0 1 1],'units','normalized');
title(['V_{bus}: ', num2str(sprintf('%0.3f', Vbus)),...
    'V;  Efficiency: ', num2str(sprintf('%4.2f', Eff)), '%']); 

drawnow;
figure(2); set(2,'position',[ss(3)/3 2*ss(4)/4+20 ss(3)/3 ss(4)/3]);   
subplot(2,1,1);
[mx,my]=meshgrid(0:nelx,0:nely);
imagesc(reshape(Vavg,nely,nelx)); axis equal; axis tight;colorbar; %shading interp;%interpolatie --> shape function;
subplot(2,1,2);
imagesc(reshape(je,nely,nelx)); axis equal; axis tight;colorbar;drawnow

% figure(3);set(3,'position',[2*ss(3)/3 2*ss(4)/4+20 ss(3)/3 ss(4)/3]); 
% subplot(2,1,1);
% [mx,my]=meshgrid(0:nelx,0:nely);
% imagesc(reshape(Pe,nely,nelx)); axis equal; axis tight;colorbar; axis([0 nelx 0 nely]);drawnow
% subplot(2,1,2);
% imagesc(reshape(Pe_loss,nely,nelx)); axis equal; axis tight;colorbar; axis([0 nelx 0 nely]);drawnow;

figure(3); set(3,'position',[6 50 ss(3)/3 ss(4)/3]);
plot(objhist(1:loop))

% figure(5); set(5,'position',[ss(3)/3 50 ss(3)/3 ss(4)/3]); 
% plot(Vbushist(1:loop))

fig=figure(1);
set(gcf,'DoubleBuffer','on'); % double buffer
set(gcf,'PaperPositionMode','auto');
M=getframe(gcf);
dummy = 0;
  end