%Initializing the parameters required by MMA optimizer
function [f0fac, f0add, m, n, xminvec, xmaxvec, low, upp, c, d, a0, a] = ...
    intialize_MMA(nelx, nely)
    f0fac = 100;          % objective scaling --> 1
    f0add = 0;          % constant addition to objective
    m = 1;                  % n of constraints / volume constraint
    n = nelx*nely;          % number of variables
    xmin = 0;               % minimum value of density (now 0)

    xminvec  = xmin*ones(n,1);  %Column vector with the lower bounds for the variables x_j.
    xmaxvec  = ones(n,1);       %Column vector with the upper bounds for the variables x_j.
    low   = xminvec;            %Column vector with the lower asymptotes from the previous  iteration (provided that iter>1).
    upp   = xmaxvec;            %Column vector with the upper asymptotes from the previous
    
    c = 1000*ones(m,1);        %Column vector with the constants c_i in the terms c_i*y_i.
    d = zeros(m,1);             %Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    a0 = 1;                     %The constants a_0 in the term a_0*z.
    a = zeros(m,1);             %Column vector with the constants a_i in the terms a_i*z
  low(n+1) = 0.5;
  upp(n+1) = 1.0;
end