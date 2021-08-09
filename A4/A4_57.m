clear;clc;
%% Q1-Q3: Explain and answer in the report
%% Q4: Fill in N1_XX.m and Ninf_XX.m in which XX is your group ID.
%% Q5: Provide your results and explanations in the report. Your code should come here to support 
H = eye(4);
lb = [2.5; -1; -2; -2]; ub = [5; 1; 2; 2];
Aeq = [1 0 -1 0;
        0 1 -0.5 -1];
Beq = [1; 0.5];
f = zeros(4,1);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[eye(4); -eye(4)],[ub; -lb],Aeq,Beq,[],[]);
%KKT conditions check
x + [eye(4); -eye(4)].'*lambda.ineqlin + Aeq.'*lambda.eqlin
% Removing lower bound
[x5l,fval5l,exitflag5l,output5l,lambda5l] = quadprog(H,f,eye(4),ub,Aeq,Beq,[],[]);

%Removing upper bound
[x5u,fval5u,exitflag5u,output5u,lambda5u] = quadprog(H,f,-eye(4),-lb,Aeq,Beq,[],[]);

% your results in the report
