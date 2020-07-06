%FUNCTION FOR DEFINING THE FREE DEGREES OF FREEDOM
function fixeddofs = fixed_dofs(nelx, nely, nroffixeddofs)


if nroffixeddofs == 1
fixeddofs = [1:nely+1];                     % left column
end
% fixeddofs = [nely-3:nely+1];                  % single point corner
% fixeddofs = [round(0.95*nely):nely+1];    %single point corner

if nroffixeddofs ==2
fixeddofs = round((nely+1)/2) + (nely+1)*round((nelx+1)/2);
end

if nroffixeddofs == 4
    fixeddofs = round((nely+1)/4) + (nely+1)*round((nelx+1)/4);
    fixeddofs = [fixeddofs fixeddofs(1)+round((nely+1)/2)];
    fixeddofs = [fixeddofs fixeddofs + (nely+1)*round((nelx+1)/2)];
end

if nroffixeddofs == 5 %one fixed busbar point at the left centre
    fixeddofs = [round((nely+1)/2)-10:round((nely+1)/2)+10];
end

if nroffixeddofs == 6
    fixeddofs = [nely-10:nely+1];
end

if nroffixeddofs == 7   % All the 4 corners
    fixeddofs = [1 nely+1 (nely+1)*(nelx+1) (nely+1)*(nelx+1)-nely];
end