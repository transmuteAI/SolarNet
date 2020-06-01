%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
function [c, dc, dv, lambda] = sens_analysis(objective, freedofs, U, K, Q, jL,...
    edofMat, KE, nelx, nely,Emin, xPhys, penal, E0, dQdv, shading,...
    penalS, Xs, Qvec, Vbus, elA, diffuse_factor)

  lambda = zeros((nely+1)*(nelx+1),1); 
  if objective == 1
  lambda(freedofs) = -K(freedofs,freedofs)\(dQdv(freedofs,freedofs)*U(freedofs)+Q(freedofs));   %(compute lambda) K = dRdU
  ce = reshape(sum((U(edofMat)*KE).*lambda(edofMat),2),nely,nelx);
  ce2 = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  cstar = (Emin+xPhys.^penal*(E0-Emin)).*ce2;
  c = Q.'*U;
%   c = U.'*Kv*U;
  dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
  if shading
      dshading = - penalS*Xs.^(penalS-1);
      dc =  -(dc + reshape(jL.*diffuse_factor.*dshading.*((U(edofMat) - lambda(edofMat))*Qvec),nely,nelx));
  end
  
  dv = ones(nely,nelx);
  
  else  
  lambda(freedofs) = -K(freedofs,freedofs)\(Vbus*sum(dQdv(freedofs,freedofs),1).');   %(compute lambda) K = dRdU
  
  ce = reshape(sum((U(edofMat)*KE).*lambda(edofMat),2),nely,nelx);
  
  c = -Vbus*sum(Q);

  dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
  if shading
      dshading = - penalS*Xs.^(penalS-1);
      dc =  -(dc + reshape(jL.*diffuse_factor.*dshading.*((-lambda(edofMat))*Qvec),nely,nelx) + Vbus*reshape(jL.*diffuse_factor.*dshading,nely,nelx)*elA);
  end
  dv = ones(nely,nelx);
   
  end