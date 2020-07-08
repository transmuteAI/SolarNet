function Auxil(nelx, nely, volfrac, penal, rmin, ft)
 %global nelx nely penal penalS rmin ft;
 penalS = 3.0; %penalty factor for shading
 disp('START');
 [Emin, E0] = cellprop(1, 1);
 jL = 310; %65; %353;       
 j0 = -0.006; %-0.7158; %-0.0029;     
 beta = 16.4; %8.038; %19.31;    %describes the p-n junction properties
 Uref = 0.7; % 0.7;%0.68;     %fixed voltage for testing
 s_Vbus = 0.7;  %factor that defines Vbus as a function of Uref
 Vbus = Uref * s_Vbus;   %voltage for fixed points in the space

 %% OPTIMIZATION OPTIONS
 shading = 1;                    %control variable for sensitivity
 jecorrect = 1; 
 Qmethod  = 3;                   %1,2 or 3
 elec = 1;                       %control variable for initial electrode
 nrofe = 991;                      %number of electrodes
 objective = 2;                  %1-> average voltage 2-> power
 optimizer = 1;                  %1-> MMA 2-> Optimality criteria
 maxiter = 1000;  
 Lx = 15e-3;                     %Length in x-direction
 elh = Lx/nelx;                  % element size
 Ly = nely*elh;                  % Length in y-direction
 elA = elh^2; 
 %disp(elA);
 ewidth = 10;
 [evec, evec2d] = create_electrode(ewidth, nrofe, nely);
 phi = compute_shape_fn();       %U*phi gives voltages....
 KE = finite_elem();
 Qvec = elA/4*ones(4,1);
 Qmat = elA/16*ones(4);
 nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
 edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
 edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);

 iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
 jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
 R = zeros((nely+1)*(nelx+1),1);
 Q = zeros((nely+1)*(nelx+1),1);
 U = zeros((nely+1)*(nelx+1),1);
 nroffixeddofs = 1;
 fixeddofs = fixed_dofs(nelx, nely, nroffixeddofs);
 nodofs = [];
 alldofs = [1:(nely+1)*(nelx+1)];
 actualdofs = setdiff(alldofs, nodofs);
 freedofs = setdiff(alldofs, fixeddofs);
 freedofs = setdiff(freedofs, nodofs);
 x = initial_design(nelx, nely, elec, 0.4, evec);
 %magesc(x);
 xPhys = x; %1 - exp(x) + x*exp(-1);
 Vavg = .25*sum(U(edofMat),2);
 %magesc(xPhys)
 Xs = 1-reshape(xPhys,nelx*nely,1);
 [je, je_check] = compute_je(Qmethod, shading, Xs, penalS, j0, jL, beta, Vavg,...
    edofMat, phi, U, jecorrect);
  %imagesc(reshape(je, 50,50));
  for i=1:size(edofMat,1)
    Q(edofMat(i,:)) = Q(edofMat(i,:))+je(i)*Qvec;
  end
  
  Edist = Emin+xPhys(:)'.^penal*(E0-Emin);
  %Edist(noelems) = 1e-6;
  sK = reshape(KE(:)*(Edist),16*nelx*nely,1);NRiter = 0;
  Qnorm = norm(Q);
  
  
  NRiterMax = 15;
  NRtol = (nelx*nely)*1e-10;
  Kv = sparse(iK,jK,sK);
  
  while 1
    NRiter = NRiter +1;
    Qnorm = norm(Q);
    if NRiter > NRiterMax, break; end    %loop control statement
    U(fixeddofs) = Vbus; %setting boundary values
    % compute residual
        R = Kv*U - Q;
    
    %norm(R(freedofs))
    if norm(R(freedofs))/Qnorm<(NRtol/100), 
%     if norm(R)/Qnorm<(NRtol/100),     
        disp(['Convergence reached in ',num2str(NRiter),' iterations']);
        break 
    else
        disp(['convergence: ' num2str(norm(R(freedofs))/Qnorm)]);
    end
    dJdV = compute_dJdV(Qmethod, j0, beta, Vavg, edofMat, phi,...
        jecorrect, je_check, U, []); 
    sKQ = reshape(Qvec(:)*dJdV(:)',16*nelx*nely,1);    
    dQdv = sparse(iK,jK,sKQ);
    K = Kv - dQdv;           %This K equals dR/dU
    
    % solve increment and update U vector;
    deltaU = K(freedofs,freedofs)\R(freedofs);
    U(freedofs) = U(freedofs) - deltaU;
    clear deltaU
    
    % compute PV currents:
    Q = zeros((nely+1)*(nelx+1),1); 
    Vavg = .25*sum(U(edofMat),2);
    
    [je, je_check] = compute_je(Qmethod, shading, Xs, penalS, j0, jL, beta, Vavg,...
    edofMat, phi, U, jecorrect);
    
    for i=1:size(edofMat,1)
    Q(edofMat(i,:))=Q(edofMat(i,:))+je(i)*Qvec;
    end  
  end 
  # Below, elA and Vbus seem ok, this means that the value of current is very high
  max_je = max(je)
  P = sum(elA*je)*Vbus              % Power of the cell = sum of all currents*voltage
  Pe = je.*Vavg*elA;    % Local Power
  Pe_eff = je*Vbus*elA; % Contribution to total power // Effective power
  Pe_loss = Pe-Pe_eff;                % Power loss
  Eff = ((P/(Lx*Ly* 1))/1000)*100;
  %disp(Lx);
  %disp(Ly);
  disp(P);
  disp(Eff);
endfunction
