  % FILTERING/MODIFICATION OF SENSITIVITIES
function [dc, dv] = filtering(ft, H, x, dc, dv, Hs, beta)
  
if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
elseif ft == 3
    dx = beta*exp(-beta*x) + exp(-beta);
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);
end
  
end