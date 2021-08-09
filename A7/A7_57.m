close all;clear;clc;
%% You can use any command in MPT in this assignment

%% Question 1
% write the code to plot X0. In the report, give Pf,  the X0, and your
% motivation for choosing this Pf.
A = [1.2 1;
    0 1];
B = [0;1]; 
Q = eye(2); R = 100;
model = LTISystem('A',A,'B',B);
model.x.min = [-15; -15];
model.x.max = [15; 15];
model.u.min = -1; model.u.max = 1;
X = Polyhedron('lb',model.x.min,'ub',model.x.max);
U = Polyhedron('lb',model.u.min,'ub',model.u.max);
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
Xf = [0;0]; 
TS = Polyhedron('Ae',eye(2),'be',Xf);
model.x.with('terminalSet');
model.x.terminalSet = TS;
N=3;
Pn = model.LQRPenalty; Pf = Pn.weight;
model.x.with('terminalPenalty');
model.x.terminalPenalty = Pn;
Xn = pre7_57(model,TS,N);
figure(1)
plot(Xn);
%% Question 2
% write the code to plot the requested figures. Provide the figures in your
% report and explain your observation about the error between the state and
% state prediction as N increases.
x0 = [4; -2.6];
N2 = 100; 
N2m = 1;
mpc = MPCController(model, N2m);
loop = ClosedLoop(mpc, model);

N= [10, 15, 20];


for i=1:numel(N)
    mpc.N = N(i);
    datasim{i} = loop.simulate(x0,N2);
    [~, ~, openloop] = mpc.evaluate(x0);
    mpceval{i} = openloop;
    
    figure('Color','White');
    hold on, grid on;
    plot(0:N2, datasim{i}.X', 'Linewidth',3);
    plot(0:N(i), mpceval{i}.X', '-.', 'Linewidth',3, 'Color','black');
    title(sprintf('N=%.f',N(i)))
    xlim([0 N2])
    xlabel 'time step', ylabel 'state position'
    legend({'x1','x2','predicted at x(0)'})

end

%% Question 3
% no code is needed. Answer in the report

%% Question 4
% write a code that calculates the figures and costs. Provide the figures
% and costs in the report. for costs, provide a table in the report that
% provides all costs for all different methods in the question (4 methods,
% each with three different costs as defined in A7 assignment). If you what
% to use some functions in the code, you can write them in different matlab
% files and submit them with the rest of your files