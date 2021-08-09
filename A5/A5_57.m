clc;clear;close all;
%% Q1
% you are supposed to find Matrix S such that V(x(k)) is a Lyapunov
% function that V is decreasing (except at the origin) in the report;
%In addition, define matrix S here as well.
A = [0.5 1; -0.1 0.2];
Sv = sym('Sv', [2,2]); Q = eye(2);
Eqn = A.'*Sv*A + Q == Sv;
Vs = solve(Eqn,Sv); 
MyFieldNames = fieldnames(Vs);
for i=1:4
Vss(i,1) = getfield(Vs,MyFieldNames{i});
end
 S= vpa([ Vss(1), Vss(2);
     Vss(3), Vss(4)]);

%% Q2
% answer in the report;
A = [1 1;
    -1 5]; 
B = [0;1]; Q = eye(2); R = 1;
[P,L,K] = dare(A,B,Q,R);

%% Q3 
% part a:
% Define N (the shortest horizon that ...) here as N3. You can use DP_XX.m
% that you have writen in previous assignments. Do note that when I run
% your code, I should be able to see the result and just writing a number
% would not be enough. Mention this N in your report as well.
Pf = eye(2);
N = 1; R=1;
K0 = zeros(1,2); 
lamb=eig(A+B*K0);
while lamb(1) > 1 || lamb(2) >1
   
    [K0,P0] = DP_57(A,B,N,Q,R,Pf);
    lamb = eig(A+B*K0);
     N = N+1;
end
N3 = N-1;


    
    


% part b:
% explain in the report



% part c:
% Fill in the function Pf_X.m in which X is your group number. Motivate
% the concept behind your code in the report.

% part d:
% you can use trial and error to find this R; just provide the R that works
% and check the stability by checking the eigenvalues of the closed loop
% system with this new R; define it in the code as Rnew
 Rnew=0.1;
 Pf = Q;
 [K03,P03] = DP_57(A,B,1,Q,Rnew,Pf);
lamb3 = abs(eig(A+B*K03));

%% Q4 
% write your proof in the report

%% Q5
% answer to the question in the report. Do note that you can verify your
% answer by checking it numerically in Matlab (this is just for your own
% and you may not provide any code regarding this)
A = [1 1; -1 5]; B = [0;1]; Q = eye(2);
Pf = Q;
syms k real
P1 = simplify(Q + A.'*Pf*A - A.'*Pf*B*inv(k + B.'*Pf*B)*B.'*Pf*A);
K5 = simplify(-inv(k + B.'*P1*B)*B.'*P1*A);

%% Q6
% answer in the report



%% do not comment this! when you run this script, the printed values should be the correct answers!
% S
% N3
% Rnew



