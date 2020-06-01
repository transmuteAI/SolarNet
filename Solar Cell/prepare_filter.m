function [H, Hs] = prepare_filter(nelx, nely, rmin)

% iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
% jH = ones(size(iH));
% sH = zeros(size(iH));
% k = 0;
% for i1 = 1:nelx
%   for j1 = 1:nely
%     e1 = (i1-1)*nely+j1;
%     for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
%       for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
%         e2 = (i2-1)*nely+j2;
%         k = k+1;
%         iH(k) = e1;
%         jH(k) = e2;
%         sH(k) = max(0.001,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
%       end
%     end
%   end
% end
% H = sparse(iH,jH,sH);
% Hs = sum(H,2);

iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = i1-(ceil(rmin)-1):i1+(ceil(rmin)-1)
            for j2 = j1-(ceil(rmin)-1):j1+(ceil(rmin)-1)
                if i2<= 0
                    e2i=0;
                elseif i2>nelx
                    e2i=(nelx-1)*nely;
                else
                    e2i=(i2-1)*nely;
                end
                if j2<=0
                    e2j=j2+nely;
                elseif j2>nely
                    e2j=j2-nely;
                else
                    e2j=j2;
                end
                e2=e2i+e2j;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0.001,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end