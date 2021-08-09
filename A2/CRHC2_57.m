function [Z,VN]=CRHC2_57(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% F1, G1, h1, F2, G2, h2 are constraint matrices
%% Be aware of the F1, F2, G1, G2, h1, and h2! 
%% x0 is the initial condition

om = A;
gam = kron(eye(N),B);

for i = 1:N-1
    om = [om; A^(i+1)];
    gam = gam + kron(diag(ones(N-i,1),-i),A^i*B);
end

Qb = blkdiag(kron(eye(N-1),Q),Pf);
Rb = kron(eye(N),R);
H = 2*(gam'*Qb*gam + Rb);
f = (2*x0'*om*Qb*gam)';
Ai = F2*gam + G2;
Bi = h2 - F2*om*x0;
Aeq = F1*gam + G1;
Beq = h1 - F1*om*xo;
end
