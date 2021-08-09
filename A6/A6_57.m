clear;clc;close all;
%% Q1
% define H and V that are correspondent to the polyhedron for H and v representation. Plot the polyhedron and
% explain the difference in the report.

A = [0 1; -1 0; -1 -1; 1 1];
b = [0; 0; 1; 1];

P1 = Polyhedron('A',A,'b',b);
P1.computeVRep

H = P1.H;
V = P1.V;


%% Q2
% define Q and P, find the sum and the difference and plot the results.
% Sketch the plots in the report as well
A1 = [0 1; 
      1 0;
      0 -1;
      -1 0];
b1 = [ 2 2 2 2 ] ;
A2 = [-1 1 1 -1;
    -1 1 -1 1].';
b2 = [1 1 1 1].';
P2 = Polyhedron('A',A1,'b',b1);
Q2 = Polyhedron('A',A2,'b',b2);
MS = plus(P2,Q2);
PD = minus(P2,Q2);
figure(1)
subplot(2,2,1);
plot(P2);
title('P');
subplot(2,2,2);
plot(Q2);
title('Q');
subplot(2,2,3);
plot(MS);
title('Minikowiski Sum');
subplot(2,2,4);
plot(PD);
title('Pontryagin Difference');

%% Q3
% write a code that shows S is invariant and explain your approach in the
% report
A = [0.8 0.4; -0.4 0.8];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P3 = Polyhedron('A',Ain,'b',bin);
figure(2)
plot(P3,'alpha',0.3);
hold on
for j = 1:1
    P3x = Polyhedron('A',Ain*inv(A)^j,'b',bin);
    plot(P3x,'alpha',0.8);
end
title('Proof showing S as positively invariant');
legend('S - Constraints', 'S - x^+ = f(x)');

%% Q4
% Fill in the Reach_X function and Plot S and its one step reachable set.
% Note that you are not supposed to change the inputs and outputs of the
% function.
A = [0.8 0.4; -0.4 0.8]; B = [0;1];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P4 = Polyhedron('A',Ain,'b',bin);
U4 = Polyhedron('lb',-1,'ub',1);
R4 = Reach_57(A,B,P4,U4);
figure(3)
plot(P4,'alpha',0.2);
hold on
plot(R4,'alpha',0.5);
legend('Set S', 'Set Reach(S)');
%% Q5
% Fill in the Pre_X function and Plot S and its Pre set.
% Note that you are not supposed to change the inputs and outputs of the
% function.
A = [0.8 0.4; -0.4 0.8]; B = [0;1];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P5 = Polyhedron('A',Ain,'b',bin);
U5 = Polyhedron('lb',-1,'ub',1);
Pr = Pre_57(A,B,P5,U5);
figure(4)
plot(P5,'alpha',0.2);
hold on
plot(Pr,'alpha',0.5);
legend('Set S', 'Set Pre(S)');
%% Q6
% part 1: Fill in the function ShorterstN_X.m and use it to find the shortest N
% that is feasible. Note that you are not supposed to change the inputs and outputs of the
% function.
A = [0.9 0.4; -0.4 0.9]; B = [0;1];
Pf = zeros(2); Q = eye(2); R = 1;
x0 = [ 2 0 ].';
xu = [3 3].'; uu = 0.1; %constraints
for N = 1:50
    [Z,exitflag]=ShortestN_57(A,B,N,Q,R,Pf,xu,uu,x0);
    if exitflag ==1
        break;
    end
end
N
% part 2: Fill in the function RHCXf_X.m and use it to check the
% feasibility. Note that you are not supposed to change the inputs and outputs of the
% function.
model = LTISystem('A',A,'B',B);
model.x.min = -xu; model.x.max = xu;
model.u.min = -uu; model.u.max = uu;
Xf2 = model.invariantSet;
N6=2;
[Z,exitflag] = RHCXf_57(A,B,N6,Q,R,Pf,xu,uu,Xf2,x0);
figure(5)
plot(Xf2,'alpha',0.3);
legend('X_{f2}');
[K0,~]=DP_57(A,B,N6,Q,R,Pf);
stability_check6=abs(eig(A+B*K0))
% part 3: Plot the feasible sets for the initial condition in part 1 and 2
% and plot those sets. Answer to the rest of the question in the report.
Xf1 = Polyhedron( 'Ae', eye(2), 'be', [0;0]);
U = Polyhedron( 'lb', -uu, 'ub', uu);

XN1 = model.reachableSet('X',Xf1,'U',U,'N',26,'direction','backward');
XN2 = model.reachableSet('X',Xf2,'U',U,'N',2,'direction','backward');

figure(6)
plot(XN2,'alpha',0.2);
hold on
plot(XN1,'alpha',0.6);

legend({'X_{N2}(N=26)','X_{N1}(N=2)'})
