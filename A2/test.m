clc;
clear all;
A=[1.0025 0.1001;0.05 1.0025];
B=[0.005;0.1001];
C=[1 0];
D=0;
Q=[5 0;0 1];
Pf=[5 0;0 1];
R=0.5;
N = 4;
[K0_5,P0_5]=BS_57(A,B,N,Q,R,Pf);
pol_5 = eig(A+B*K0_5(1));

N5=4;