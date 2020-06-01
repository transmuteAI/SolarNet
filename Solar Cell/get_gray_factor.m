function gray_factor = get_gray_factor(xPhys1, nelx, nely, x_low)

vol_frac = mean(xPhys1(:));

xPhys2 = xPhys1 > x_low;
elem_frac = sum(xPhys2(:))/(nely*nelx);

gray_factor = 1 - (vol_frac/elem_frac);
clear vol_frac xPhys1 xPhys2;

end