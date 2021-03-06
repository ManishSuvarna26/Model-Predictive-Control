function [x,fval]=N1_57(A,b)
%% This function finds x to minimize the first norm
%% The inputs are Matrix a and vector b
%% The outputs are the optimal x and the norm value (fval)
%% Do not change the inputs and outputs!
%% use linprog to solve the problem
n = length(b);
f = [zeros(n,1); ones(n,1)];
Ai = [A -eye(n);
      -A -eye(n)];
bi = [b; -b];
[x,fval] = linprog(f,[],[],Ai,bi);


end

