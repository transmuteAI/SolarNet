%FUNCTION FOR COMPUTATION OF SHAPE FUNCTIONS
function phi = compute_shape_fn()
Xa = [0.211324865405187 0.788675134594813 0.788675134594813 0.211324865405187];
Yb = [0.211324865405187 0.211324865405187 0.788675134594813 0.788675134594813];
a = 1 ; b = 1;

phi1 = (1-Xa./a).*(Yb./b);
phi2 = (1-Xa./a).*(1-Yb./b);
phi3 = (Xa./a).*(Yb./b);
phi4 = (Xa./a).*(1-Yb./b);

%put them in the order 2 4 3 1!!
phi = [phi2;phi4;phi3;phi1];
end
