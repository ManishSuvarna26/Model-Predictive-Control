function [K0,P0]=DP_57(A,B,N,Q,R,Pf)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% K0 is the controller gain when u(0)=K0x(0)
%% P0 describes the final cost as VN=(x0^T)*P0*x0 
P0 = Pf;
for k = N:-1:1
    K0 = -inv(R + B'*P0*B)*B'*P0*A;
    P0 = Q + A'*P0*A - A'*P0*B*inv(R+B'*P0*B)*B'*P0*A;
end
end
