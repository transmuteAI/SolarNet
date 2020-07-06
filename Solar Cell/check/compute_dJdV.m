%FUNCTION TO COMPUTE dJ/dV
function dJdV = compute_dJdV(Qmethod, j0, beta, Vavg, edofMat, phi,...
    jecorrect, je_check, U, dofs)
 
if Qmethod == 1
    dJdV = repmat(1/4*(j0*beta*exp(beta*Vavg(:)')),4,1);
elseif Qmethod == 2
    dJdV = 1/4*(j0*beta*exp(beta*U(edofMat)))';
else
    %U(dofs) = 0;
    dJdV = 1/4*phi*(j0*beta*exp(beta*U(edofMat)*phi))';
end

if jecorrect == 1
dJdV(:,je_check) = 0;                                       %makes derivative zero if je < 0  trick!S       
end