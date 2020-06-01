function [evec, evec2d] = create_electrode(ewidth, nrofe, nely)
%function for creating electrode
if nrofe == 0
    evec = 0;
    evec2d = 0;
elseif nrofe == 1
    elnr = nely/(nrofe+1);
    evec = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr;
    evec2d = 0;
elseif nrofe == 2;
    elnr = nely/(nrofe*2);
    evec = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr;
    elnr2 = 3*nely/(nrofe*2);
    evec2 = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr2;
    evec = [evec evec2];
    evec2d = 0;
elseif nrofe == 991
    evec = 1:2:nely;
    %evec = [evec evec-1 evec+1];
    %evec = [evec-2 evec-1 evec evec+1 evec+2];
    evec2d = 0;
elseif nrofe == 10001
    evec = 30:859;
    evec2d = 210:225;
elseif nrofe == 881
    evec = floor(nely * rand(1));
    evec2d = 0;
    
elseif nrofe == 12121
    evec = 299:301;
    evec2d = 0;
end
end