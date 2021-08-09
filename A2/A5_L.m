clc;clear;close all;
%% Q1
% you are supposed to find Matrix S such that V(x(k)) is a Lyapunov
% function that V is decreasing (except at the origin) in the report;
%In addition, define matrix S here as well.

A = [0.5 1; -0.1 0.2];
Q = eye(2);

S= dlyap(A,Q)

%% Q2

A = [1 1; -1 5];
B = [0; 1];
Q = eye(2);
R = 1;

[P,lambda,K] = dare(A,B,Q,R)



%% Q3 
% part a:
% Define N (the shortest horizon that ...) here as N3. You can use DP_XX.m
% that you have writen in previous assignments. Do note that when I run
% your code, I should be able to see the result and just writing a number
% would not be enough. mention this N in your report as well.

Pf = Q;

stable = false;
N3 = 0;
while ~stable
    N3 = N3+1;
    [K0,P0] = DP_09(A,B,N3,Q,R,Pf);
    if all(abs(eig(A-B*K0))<=1)
        stable = true;
    end
end

Pf
N3
eig(A-B*K0)
P0
K0




% part b:
N = 1;
for i=1:10
    [K0,P0] = DP_09(A,B,N,Q*i^3,R,Pf);
    P0
    Q*i^3
    K0
    abs(eig(A-B*K0))
end


% part c:
% Fill in the function Pf_XX.m in which XX is your group number. Motivate
% the concept behind your code in the report.
Pf = Pf_09(A,B,Q,R)
N = 1;
[K0,P0] = DP_09(A,B,N,Q,R,Pf);
abs(eig(A-B*K0))


% part d:
% you can use trial and error to find this R; just provide the R that works
% and check the stability by checking the eigenvalues of the closed loop
% system with this new R; define it in the code as Rnew

Pf = Q;
N = 1;
Rnew = 1;

stable = false;
while ~stable && ~isinf(Rnew)
    [K0,P0] = DP_09(A,B,N,Q,Rnew,Pf);
    if all(abs(eig(A-B*K0))<=1)
        stable = true;
        break;
    end
    Rnew = Rnew/10
end

Pf
Rnew
eig(A-B*K0)
P0
K0


%% Q4 
% write your proof in the report
%% Q5
% answer to the question in the report. Do note that you can verify your
% answer by checking it numerically in Matlab (this is just for your own
% and you should not provide any code regarding this)

% PART A

syms Rs
K = @(R) (R + B'*Pf*B) \ B'*Pf*A;
eigT = @(R) eig(A-B*K(R))

eigT(Rs)

abs(eigT(0.091))

Rv = linspace(0.001,2,200);
for i=1:numel(Rv)
    polesabs(i,:) = abs(eigT(Rv(i)));
end
figure('Color','white'), grid on, hold on;
plot(Rv, polesabs, 'Linewidth', 2)
xlabel 'R', ylabel 'abs(eig(A-BK))', title 'N=1'


% PART B

P1 = @(R) Q + A'*Pf*A - A'*Pf*B*(R + B'*Pf*B)^-1*B'*Pf*A;
P1(Rs)
K = @(R) (R + B'*P1(R)*B) \ B'*P1(R)*A;
eigT = @(R) eig(A-B*K(R))

simplify(eigT(Rs))

Rv = linspace(0.001,10,200);
for i=1:numel(Rv)
    polesabs(i,:) = abs(eigT(Rv(i)));
end
figure('Color','white'), grid on, hold on;
plot(Rv, polesabs, 'Linewidth', 2)
xlabel 'R', ylabel 'abs(eig(A-BK))', title 'N=2'

%% Q6
% answer in the report