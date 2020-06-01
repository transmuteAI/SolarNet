function dummy = save_output(nelx, nely, U, Q, xPhys, j0, jV, P, Eff, beta, Vbus, elA, c, filename)

 jemax = j0+jV*(exp(beta*Vbus)-1);   % Current density if cell has zero resistance
 Pmax = elA*jemax*(nelx*nely)*Vbus   % Maximum power possible (withouth shading etc)
 Umax = max(U);                      % maximum voltage found
 je_Umax = j0+jV*(exp(beta*Umax)-1);
  
 disp(['Maximum possible power from cell: Pmax: ' num2str(Pmax)])
 disp(['Fraction obtained: P/Pmax :  ' num2str(P/Pmax)])
 disp(['Fraction obtained: Pshading/Pmax :  ' num2str(c/Pmax)])
 disp(['Maximum Voltage found: ' num2str(max(U))])
 disp(['Corresponding current: ' num2str(je_Umax)])
 
 save([filename, 'workspace']);
 Q = reshape(Q, (nely+1), (nelx+1));
 U = reshape(U, (nely+1), (nelx+1));
 xPhys = reshape(xPhys, nely, nelx);
 save([filename, 'current.dat'], 'Q', '-ascii');
 save([filename, 'voltage.dat'], 'U', '-ascii');
 save([filename, 'density.dat'], 'xPhys', '-ascii');
 
end