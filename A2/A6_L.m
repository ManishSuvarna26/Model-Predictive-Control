clear;clc;close all;
%% Q1
% define H and V that are correspondent to the polyhedron for H and v representation. Plot the polyhedron and
% explain the difference in the report.

% mpt_demo_sets1
% mpt_demo2

A = [0 1; -1 0; -1 -1; 1 1];
b = [0; 0; 1; 1];

P1 = Polyhedron('A',A,'b',b);
P1.computeVRep

H = P1.H
V = P1.V

% plot(P1)

% latex(sym(H))
% latex(sym(V))


%% Q2
% define Q and P, find the sum and the difference and plot the results.
% Sketch the plots in the report as well

A1 = [0 1; 1 0; 0 -1; -1 0];
b1 = [2;2;2;2];
A2 = [-1 -1; 1 1; 1 -1; -1 1];
b2 = [1;1;1;1];

P = Polyhedron('A',A1,'b',b1)
Q = Polyhedron('A',A2,'b',b2)

PpQ = plus(P,Q)
PmQ = minus(P,Q)

figure('Color','white','Position',[387   356   744   126])
subplot(1,4,1)
plot(P)
title 'P'
subplot(1,4,2)
plot(Q)
title 'Q'
subplot(1,4,3)
plot(PpQ)
title 'P+Q'
subplot(1,4,4)
plot(PmQ)
title 'P-Q'
% fp.savefig('mink')

% figure
% plot(P,'alpha',0.2,Q,'alpha',0.2,PpQ,'alpha',0.2,PmQ,'alpha',0.2)


%% Q3
% write a code that shows S is invariant and explain your approach in the
% report

A = [0.8 0.4; -0.4 0.8];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P = Polyhedron('A',Ain,'b',bin);


figure('Color','white','Position',[573   416   248   227])
plot(P,'alpha',0.2);
hold on;
for i=1:1
    P_R = Polyhedron('A',Ain*inv(A)^i,'b',bin);
    plot(P_R,'alpha',0.6);
end
axis equal
legend({'S','Reach(S)'})
% fp.savefig('cis')



%% Q4
% Fill in the Reach_XX function and Plot S and its one step reachable set.
% Note that you are not supposed to change the inputs and outputs of the
% function.

A = [0.8 0.4; -0.4 0.8];
B = [0;1];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P = Polyhedron('A',Ain,'b',bin);

U = Polyhedron('lb',-1,'ub',1);

Reach = Reach_09(A,B,P,U);

figure('Color','white','Position',[608   390   179   336])
plot(P,'alpha',0.6,Reach,'alpha',0.2)
axis equal
legend({'S','Reach(S)'})
% fp.savefig('reach')



%% Q5
% Fill in the Pre_XX function and Plot S and its Pre set.
% Note that you are not supposed to change the inputs and outputs of the
% function.

A = [0.8 0.4; -0.4 0.8];
B = [0;1];

Ain = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
P = Polyhedron('A',Ain,'b',bin);

U = Polyhedron('lb',-1,'ub',1);

Pre = Pre_09(A,B,P,U);

figure('Color','white','Position',[603   514   186   207])
plot(P,'alpha',0.6,Pre,'alpha',0.2)
axis equal
legend({'S','Pre(S)'})
% fp.savefig('pre')


%% Q6
% part 1: Fill in the function ShorterstN_XX.m and use it to find the shortest N
% that is feasible. Note that you are not supposed to change the inputs and outputs of the
% function.

A = [0.9 0.4 ; -0.4 0.9];
B = [0;1];
Pf = zeros(2);
Q = eye(2);
R=1;
x0 = [2;0];

x_ub = [3;3];
u_ub = 0.1;


for N=1:100
    [Z,exitflag] = ShortestN_09(A,B,N,Q,R,Pf,x_ub,u_ub,x0);
    if exitflag==1
        break;
    end
end
N



% part 2: Fill in the function RHCXf_XX.m and use it to check the
% feasibility. Note that you are not supposed to change the inputs and outputs of the
% function.

model = LTISystem('A', A, 'B', B);
model.x.min = -x_ub;
model.x.max = x_ub;
model.u.min = -u_ub;
model.u.max = u_ub;
Xf2 = model.invariantSet;
% plot(Xf)

N=2;
[Z,exitflag] = RHCXf_09(A,B,N,Q,R,Pf,x_ub,u_ub,Xf2,x0)

figure('Color','white', 'Position',[454   448   311   298])
plot(Xf2,'color','red','alpha',0.5)
axis equal
legend({'X_f2'})
% fp.savefig('xf2')


% Check if system is stable for x0 \in Xf
[K0,P0]=DP_09(A,B,N,Q,R,Pf)
abs(eig(A-B*K0))



% part 3: Plot the feasible sets for the initial condition in part 1 and 2
% and plot those sets. Answer to the rest of the question in the report. 


Xf1 = Polyhedron( 'Ae', eye(2), 'be', [0;0]);
U = Polyhedron( 'lb', -u_ub, 'ub', u_ub);

XN1 = model.reachableSet('X',Xf1,'U',U,'N',26,'direction','backward')
XN2 = model.reachableSet('X',Xf2,'U',U,'N',2,'direction','backward')

figure('Color','white', 'Position',[454   448   311   298])
plot(XN2,'color','red','alpha',0.2, XN1,'color','red','alpha',0.6)
axis equal
legend({'X_N2','X_N1'})
% fp.savefig('xn')