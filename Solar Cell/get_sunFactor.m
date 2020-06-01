% This function designs the illumination pattern on the solar cell

function sunFactor = get_sunFactor(pType, nely, nelx)

    sunFactor = ones(nely, nelx);
    
    if pType == 1
        sunFactor(:, fix(nelx/2)+1:nelx) = 50;
    end
    
    imagesc(sunFactor);
end