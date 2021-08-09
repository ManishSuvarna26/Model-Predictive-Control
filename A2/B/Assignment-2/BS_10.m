function [K0,P0]=BS_10(A,B,N,Q,R,Pf)
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% K0 is the controller gain when u(0)=K0x
%% P0 describes the final cost as VN=(x0^T)*P0*x0 
P0=Pf
Omg = A;
for i=2:N
    Omg =[Omg;A^i];
end
Gam =zeros(N*2,N);
for i=0:N-1
    T_v =(A^i)*B;   
    Gam = kron(diag(ones(N-i,1),-i),T_v)+Gam;
end
Q_diag= blkdiag(Q,Q,Q,Pf);
R_diag=blkdiag(R,R,R,R);
    K0 = -inv((Gam'*Q_diag*Gam)+R_diag)*(Gam'*Q_diag*Omg);
    P0 = (Q+(Omg'*Q_diag*Omg)-(Omg'*Q_diag*Gam)*inv((Gam'*Q_diag*Gam)+R_diag)*(Gam'*Q_diag*Omg));
end