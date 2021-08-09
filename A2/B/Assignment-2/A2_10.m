.clear;
clc;
%% Q1: Fill the DP_XX.m function using the dynamic programic approach.

%% 2: Find the shortest N that stabilizes the system using DP_XX.m ; define N2 that gives the shortest N, i.e. N2=min(N) subject to the system stability
A=[1.0025 0.1001;0.05 1.0025];  
B=[0.005;0.1001];
C=[1 0];
D=0;
Q=[5 0;0 1];
Pf=[5 0;0 1];
R=0.5;
N=4;
[K0,P0]=DP_10(A,B,N,Q,R,Pf);
N2=4;
Eigen_0=eig(A+B*K0);

%% Q3: Define Pf3 as the solution to the Riccati Equation; define N3 that gives the shortest N, i.e. N3=min(N) subject to the system stability
Pf3 = dare(A,B,Q,R);
Pf=Pf3;
[K0,P0]=DP_10(A,B,N,Q,R,Pf);
N3=1;
Eigen_1=eig(A+B*K0);

%% Q4: Fill the BS_XX.m function using the batch solution approach.
[K0,P0]=BS_10(A,B,N,Q,R,Pf);

%% Q5: Find the shortest N that stabilizes the system using BS_XX.m; define N5 that gives the shortest N, i.e. N5=min(N) subject to the system stability

Eigen_2=eig(A+B*K0(1));
N5=4;

%% Q6: Use BS_XX.m or DP_XX.m and simulate the system for 20 steps; plot the inputs and the states for these four cases.
Pf=[5 0;0 1];
tf=20;
x0=[1;0];
R=0.5;
N=5;
[K0]=DP_10(A,B,N,Q,R,Pf);
t= 0:0.1:tf;
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
f1=plot(t,Y_1,'k');
hold on 
figure(2)
f2=plot(t,U_1,'k');
hold on

R=0.5;N=15;
[K0]=DP_10(A,B,N,Q,R,Pf);
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
f3=plot(t,Y_1,'g');
hold on
figure(2)
f4=plot(t,U_1,'g');
hold on

R=0.05;N=5;
[K0]=DP_10(A,B,N,Q,R,Pf);
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
f5=plot(t,Y_1,'r');
hold on 
figure(2)
f6=plot(t,U_1,'r');
hold on

R=0.05;N=15;
[K0]=DP_10(A,B,N,Q,R,Pf);
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
title('System Output Vs Time');
xlabel('Time in S');
ylabel('System Output');
f7=plot(t,Y_1,'m');
hold on 
figure(2);
title('System Input Vs Time');
xlabel('Time in S');
ylabel('System Input');
f8=plot(t,U_1,'m');
P1=[f1(1),f3(1),f5(1),f7(1)];
P2=[f2(1),f4(1),f6(1),f8(1)];
legend(P1,'R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15');
legend(P2,'R=0.5 & N=5','R=0.5 & N=15','R=0.05 & N=5','R=0.05 & N=15');
%% Q7: Fill the CRHC1_XX.m function
%% Q8: Fill the CRHC2_XX.m function
%% Q9: Solve Q6 using CRHC1_XX.m or CRHC2_XX.m considering the given constraints for 100 sample times
tf=100;
x0=[1;0];
% R=0.5;N=5;
% R=0.5;N=15;
% R=0.05;N=5;
% R=0.05;N=15;



