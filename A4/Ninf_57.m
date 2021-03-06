function [x,fval]=Ninf_XX(A,b)
%% This function finds x to minimize the infinity norm
%% The inputs are Matrix a and vector b
%% The outputs are the optimal x and the norm value (fval)
%% Do not change the inputs and outputs!
%% use linprog to solve the problem
n = length(b);
f = [zeros(n,1); ones(n,1)];
Ai = [A -ones(n);
      -A -ones(n)];
bi = [b; -b];
[x,fval] = linprog(f,[],[],Ai,bi);

end