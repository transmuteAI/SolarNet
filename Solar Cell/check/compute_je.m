%FUNCTION TO COMPUTE ELEMENT CURRENT DENSITY
function [je, je_check] = compute_je(Qmethod, shading, Xs, penalS, j0, jL, beta,...
    Vavg, edofMat, phi, U, jecorrect)

if Qmethod == 1
    if shading
    je = Xs.^penalS*jL+j0*(exp(beta*Vavg)-1);                 %current density for 1 node
    else
    je = jL+j0*(exp(beta*Vavg)-1);
    end
elseif Qmethod == 2
    if shading
    je = Xs.^penalS*jL+1/4*sum(j0*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    else
    je = jL+1/4*sum(j0*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    end
else
    if shading
    je = Xs.^penalS*jL+1/4*sum(j0*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    else
    je = jL+1/4*sum(j0*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    end
end


disp(['Highest voltage found during iteration: ' num2str(max(U)) ]) 
disp(['Lowest current found during iteration: ' num2str(min(je)) ]) 
je_check = je < (-jL);

if any(je_check)
    disp('Voltage above upper limit, solution not correct')
    if jecorrect == 1
        je(je_check) = 0;
        disp('Correction on current density applied')
    end
end
end
