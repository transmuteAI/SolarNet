%function for making Vbus as a parameter
function dc_Vbus = sens_analysis_Vbus(lambda, Qmethod, j0, beta, Vavg, edofMat,...
    phi,jecorrect, je_check, U, Uref, Q, s_Vbus, fixeddofs, freedofs,...
    Kv, Qvec, nelx, nely, iK, jK)
  dJfdV = compute_dJdV(Qmethod, j0, beta, Vavg, edofMat, phi,...
         jecorrect, je_check, U, fixeddofs); 
  sKQ1 = reshape(Qvec(:)*dJfdV(:)',16*nelx*nely,1);
  dQfdvf = sparse(iK,jK,sKQ1);
  dQfdvf = dQfdvf(freedofs, freedofs);
  dQfdvp = sparse(iK, jK, sKQ1);
  dQfdvp = dQfdvp(freedofs, fixeddofs);
  
  %computing dQdv
  dJdV1 = compute_dJdV(Qmethod, j0, beta, Vavg, edofMat, phi,...
         jecorrect, je_check, U, []); 
  sKQ = reshape(Qvec(:)*dJdV1(:)',16*nelx*nely,1);
  dQdv1 = sparse(iK,jK,sKQ);
  
  dQdvf = dQdv1(freedofs, freedofs);
  dfdVf = s_Vbus*Uref*sum(dQdvf);
  dQdvp = dQdv1(fixeddofs, fixeddofs);
  dfdVp = s_Vbus*Uref*sum(dQdvp);
  %lambdaf = -(Kv(freedofs, freedofs) - dQfdvf)\dfdVf';
  lambdaf = lambda(freedofs);
  dfds_Vbus = Uref*sum(Q);
  dUpds_Vbus = Uref;
  
  dctemp = dfdVp + lambdaf'*Kv(freedofs, fixeddofs) - lambdaf'*dQfdvp;
  dc_Vbus = dfds_Vbus + sum(dctemp*dUpds_Vbus);
  dc_Vbus = dc_Vbus*-1;
end