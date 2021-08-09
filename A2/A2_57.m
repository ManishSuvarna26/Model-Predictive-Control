clear;clc;
%% Note that when you run this file, N2, N3, Pf3, and N5 should be defined in Workspace and have the correct values.
%% Q1: Fill the DP_57.m function using the dynamic programic approach.

%% Q2: Find the shortest N that stabilizes the system using DP_57.m ; define N2 that gives the shortest N, i.e. N2=min(N) subject to the system stability
A=[1.0025 0.1001;0.05 1.0025];
B=[0.005;0.1001];
C=[1 0];
D=0;
Q=[5 0;0 1];
Pf=[5 0;0 1];
R=0.5;
N=4;
[K0,P0] = DP_57(A,B,N,Q,R,Pf);
pol_0 = eig(A+B*K0);

N2= 4;

%% Q3: Define Pf3 as the solution to the Riccati Equation; define N3 that gives the shortest N, i.e. N3=min(N) subject to the system stability
[X,L,G] = dare(A,B,Q,R);
Pf3 = X;
N = 1;
[K0_3,P0_3]=DP_57(A,B,N,Q,R,Pf3);
pol_3 = eig(A+B*K0_3);

N3 = 1;
%Pf3= the shortest N!

%% Q4: Fill the BS_XX.m function using the batch solution approach.
N = 4;
[K0_5,P0_5]=BS_57(A,B,N,Q,R,Pf);

%% Q5: Find the shortest N that stabilizes the system using BS_XX.m; define N5 that gives the shortest N, i.e. N5=min(N) subject to the system stability


pol_5 = eig(A+B*K0_5(1));

N5=4;

%% Q6: Use BS_XX.m or DP_XX.m and simulate the system for 20 steps; plot the inputs and the states for these four cases.
tf=20;
x0=[1;0];

%R=0.5;N=5;
% R=0.5;N=15;
% R=0.05;N=5;
% R=0.05;N=15;
R=0.5;
N=5;
[K0]=DP_57(A,B,N,Q,R,Pf);
t= 0:1:tf;
X_1= zeros(length(t),2);
X_1(1,:)= x0;
Y_1= zeros(length(t));
U_1= zeros(length(t));
for i=1:length(t)   
    U_1(i)=K0*X_1(i,:)';
    Y_1(i)=C*X_1(i,:)';
    X_1(i+1,:)=A*X_1(i,:)'+B*U_1(i)';
end
figure(1)
f1=plot(t,Y_1,'k','LineWidth',1);
hold on 
figure(2)
f2=plot(t,U_1,'k','LineWidth',1);
hold on

R=0.5;N=15;
[K0]=DP_57(A,B,N,Q,R,Pf);
t= 0:1:tf;
x0=[1;0];
X_1= zeros(length(t),2);
X_1(1,:)= x0;
Y_1= zeros(length(t));
U_1= zeros(length(t));
for i=1:length(t)
    U_1(i)=K0*X_1(i,:)';
    Y_1(i)=C*X_1(i,:)';
    X_1(i+1,:)=A*X_1(i,:)'+B*U_1(i)';
end
figure(1)
f3=plot(t,Y_1,'g','LineWidth',1);
hold on
figure(2)
f4=plot(t,U_1,'g','LineWidth',1);
hold on

R=0.05;N=5;
[K0]=DP_57(A,B,N,Q,R,Pf);
t= 0:1:tf;
x0=[1;0];
X_1= zeros(length(t),2);
X_1(1,:)= x0;
Y_1= zeros(length(t));
U_1= zeros(length(t));
for i=1:length(t)
    U_1(i)=K0*X_1(i,:)';
    Y_1(i)=C*X_1(i,:)';
    X_1(i+1,:)=A*X_1(i,:)'+B*U_1(i)';
end
figure(1)
f5=plot(t,Y_1,'r','LineWidth',1);
hold on 
figure(2)
f6=plot(t,U_1,'r');
hold on

R=0.05;N=15;
[K0]=DP_57(A,B,N,Q,R,Pf);
t= 0:1:tf;
x0=[1;0];
X_1= zeros(length(t),2);
X_1(1,:)= x0;
Y_1= zeros(length(t));
U_1= zeros(length(t));
for i=1:length(t)
    U_1(i)=K0*X_1(i,:)';
    Y_1(i)=C*X_1(i,:)';
    X_1(i+1,:)=A*X_1(i,:)'+B*U_1(i)';
end
figure(1);
title('System Output Vs Time(RHC)');
xlabel('Time in S');
ylabel('y(k)');
f7=plot(t,Y_1,'m','LineWidth',1);
grid on;
hold on 
figure(2);
title('System Input Vs Time(RHC)');
xlabel('Time in S');
ylabel('u(k)');
grid on;
f8=plot(t,U_1,'m','LineWidth',1);
P1=[f1(1),f3(1),f5(1),f7(1)];
P2=[f2(1),f4(1),f6(1),f8(1)];
legend(P1,'R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15','Location','northeast');
legend(P2,'R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15','Location','southeast');

%% Q7: Fill the CRHC1_X.m function
%% Q8: Fill the CRHC2_X.m function
%% Q9: Solve Q6 using CRHC1_X.m or CRHC2_X.m considering the given constraints for 100 sample times
tf=100;
x0=[1;0];
% R=0.5;N=5;
% R=0.5;N=15;
% R=0.05;N=5;
% R=0.05;N=15;
R = [0.5, 0.5, 0.05, 0.05];
N = [5, 15, 5, 15];
F1N = zeros(2);
G1N = [0;0];
h1N = zeros(2,1);
F2N = [0 1;
       0 -1;
       0 0;
       0 0];
G2N = [0;
        0;
        1;
       -1];
h2N = [0.5;
        0.5;
        0.7;
        0.7];
Q91 = figure('Color','white'); hold on, grid on;
xlabel 'time-step k', ylabel 'x(k)'
Q92 = figure('Color','white'); hold on, grid on;
xlabel 'time-step k', ylabel 'u(k)'
clr = lines(20);
for i=1:numel(R)
    F2 = kron(eye(N(i)),F2N);
    G2 = kron(eye(N(i)),G2N);
    h2 = kron(ones(N(i),1),h2N);
    F1 = kron(eye(N(i)),F1N);
    G1 = kron(eye(N(i)),G1N);
    h1 = kron(ones(N(i),1),h1N);

    [x91,u91] = Q9(A,B,N(i),Q,R(i),Pf,F1,G1,h1,F2,G2,h2,x0,tf);
    figure(Q91)
    plot(0:tf,x91(1,:)','Color',clr(i,:),'LineWidth',1,'DisplayName',sprintf('x(1) R=%.2f N=%d',R(i),N(i)));
    plot(0:tf,x91(2,:)','Color',clr(i,:),'LineWidth',1,'DisplayName',sprintf('x(2)'));
    figure(Q92)
    plot(0:(tf-1),u91','Color',clr(i,:),'LineWidth',1,'DisplayName',sprintf('u R=%.2f N=%d',R(i),N(i)));
    
end
   
figure(Q91)
legend('R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15','Location','southeast');
title('System output vs Time (CRHC)');

figure(Q92)
legend('R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15','Location','southeast');
title('System input vs Time (CRHC)');
%legend('Location','southeast');
% fp.savefig('q9_u');
function [x,u] = Q9(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0,tf)
    x = x0;
    u = [];
    for t=1:tf
        [Z,~]=CRHC1_57(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x(:,t));
        u(:,t) = Z(end-N+1:end-N+size(B,2),:);
        x(:,t+1) = A*x(:,t) + B*u(:,t);
    end
end
