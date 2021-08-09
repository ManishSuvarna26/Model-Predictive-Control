function [K0,P0]=BS_57(A,B,N,Q,R,Pf)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% K0 is the controller gain when u(0)=K0x
%% P0 describes the final cost as VN=(x0^T)*P0*x0 
%   P0=Pf
om = A;
gam = kron(eye(N),B);

for i = 1:N-1
    om = [om; A^(i+1)];
    gam = gam + kron(diag(ones(N-i,1),-i),A^i*B);
end

Qb = blkdiag(kron(eye(N-1),Q),Pf);
Rb = kron(eye(N),R);
K0 = -inv(Rb + gam'*Qb*gam)*gam'*Qb*om;
P0 = Q + om'*Qb*om - om'*Qb*gam*inv(Rb + gam'*Qb*gam)*gam'*Qb*om;