function [Emin, E0] = cellprop(cellid, subid)

% cellid = 1 : thin film Si solar cell
% cellid = 2 : organic solar cell

subid = 1;

if cellid == 1
    %TCO layer
    sigmaTL = 1e5;
    hTL = 200e-9;
    %electrode layer
    sigmaEL = 1e7;
    hEL = 10e-6;
    
elseif cellid == 2
    if subid == 1
        sigmaTL = 5e4;
        hTL = 2e-7;
        
        sigmaEL = 6.3e6;
        hEL = 3e-7;
    end
end

Emin = sigmaTL * hTL;
E0 = sigmaEL * hEL;

Emin
E0