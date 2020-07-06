%Prepare Finite Element Analysis
function KE = finite_elem()
A = [ 2/3 -1/6
     -1/6  2/3];
B = [-1/3 -1/6 
     -1/6 -1/3]; 
KE = [A B;B A];
end