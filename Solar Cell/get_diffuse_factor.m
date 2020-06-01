function diffuse_factor = get_diffuse_factor(diffuse_factor, nely, nelx, design)
xc = nelx;
yc = nely;

dmax = sqrt(((nelx-1)^2 + (nely-1)^2))*0.5;

i = 1:1:length(diffuse_factor)
xi = floor(i/(nely));
yi = mod(i, nely);
diffuse_factor = sqrt(abs(xi - xc).^2 + abs(yi - yc).^2);
diffuse_factor =  1 + ( 40*(diffuse_factor ./ dmax));

%diffuse_factor(nelx/2 * nely +1 : nelx*nely) = 1;
%diffuse_factor(1: nelx/2 * nely) = 1;

end
