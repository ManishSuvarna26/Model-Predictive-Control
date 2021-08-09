clear;
close all;
clc;
%% Question 1
A=diag([0.5 0.6 0.5 0.6]);
B=[diag([0.5 0.4]);diag([0.25 0.6])];
C=[1 1 0 0;0 0 1 1];
z_sp=[1;-1];

% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report
A1 = [eye(4)-A, -B; C, zeros(2)];
sp = [zeros(4,1); z_sp];
spt = inv(A1)*sp;
xs1= spt(1:4,1); 
us1= spt(5:6,1);
% Regulation check:
ysc = C*xs1;
%% Question 2
A=diag([0.5 0.6 0.5 0.6]);
B=[0.5;0;0.25;0];
C=[1 1 0 0;0 0 1 1];
z_sp=[1;-1];
% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report
Aeq = [eye(4) - A, -B];
Beq = zeros(4,1);
Q = eye(2);
H = 2*blkdiag(C'*Q*C,zeros(1));
f = [(-2*z_sp'*Q*C)'; zeros(1)];
x2s = quadprog(H,f,[],[],Aeq,Beq);

xs2= x2s(1:4);
us2= x2s(5);
%% Question 3
A=diag([0.5 0.6 0.5 0.6]);
B=[diag([0.5 0.4]);diag([0.25 0.6])];
C=[1 1 0 0];
z_sp=1;

% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report
H = 2*blkdiag(zeros(4),eye(2));
f = zeros(6,1);
Aeq3 = [eye(4)-A, -B;
        C, zeros(1,2)];
Beq3 = [zeros(4,1);z_sp];
x3s = quadprog(H,f,[],[],Aeq3,Beq3);
xs3= x3s(1:4);
us3= x3s(5:6);
%% second part
tf=50;                
%==========================================================================
% Process model
%==========================================================================

h = 1; %sampling time in minutes

A = [ 0.2681   -0.00338   -0.00728;
      9.703    0.3279   -25.44;
         0         0       1   ];
B = [ -0.00537  0.1655;
       1.297   97.91 ;
       0       -6.637];
C = [ 1 0 0;
      0 1 0;
      0 0 1];
Bp = [-0.1175;
      69.74;
       6.637 ];
   
n = size(A,1); % n is the dimension of the state
m = size(B,2); % m is the dimension of the control signal
p = size(C,1); % p is the dimension of the measured output

d=0.01*[zeros(1*tf/5,1);ones(4*tf/5,1)]; %unmeasured disturbance trajectory

x0 = [0.01;1;0.1]; % initial condition of system's state

%==========================================================================
% Observer model
%==========================================================================
%% choose one case, i.e. a, b, or c and then write the code for that case! for the other ones
% you just need to change the example case!
example = 'c';
switch example
    case 'a'
        nd = 2;
        Bd = zeros(n,nd);
        Cd = [1 0;0 0; 0 1]; 
    case 'b'
        nd=3;
        Bd = zeros(n,nd); 
        Cd = [1 0 0;0 0 1;0 1 0];
    case 'c'
        nd=3;
        %Bd = [0 0 0.1655;0 0 97.91; 0 0 -6.637]; 
        Bd = [zeros(3,2) Bp];
        Cd = [1 0 0;0 0 0;0 1 0];
end

%% Question 4
%Augment the model with constant disturbances; check the detectability for case "a", "b", and "c" and
%report the detectability of each one in the report

% define Ae, Be, and Ce which are the matrices for the augmented system. No
% need to report these in the report

Ae = [A, Bd;
      zeros(nd,n), eye(nd);];

Be = [B;zeros(nd,m)];

Ce = [C Cd];
% Detectability condition check:
rank_aug = rank([eye(n) - A, -Bd;
                C, Cd;])

%% Question 5
% Calculate Kalman filter gain and name it Le; no need to report this in
% the report
Q = eye(n+nd);
R = eye(p);
[X,L,G] = dare(Ae',Ce',Q,R);
Le = X*Ce'*inv(Ce*X*Ce' + R);
if example == 'a'
    Le_a = Le;
end
if example == 'c'
    Le_c = Le;
end
%% Question 6
%Target Selector
H = [1 0 0;0 0 1]; 

% Matrices for steady stae target calculation

Mss = inv([eye(n) - A, -B; H*C, zeros(size(H*C,1),m)])*[Bd; -H*Cd];
%note that Mss is defined by [xs;us]=Mss*d and you will use it later on; no need to report this in
% the report

%% Question 7
%==========================================================================
% Setup MPC controller
%==========================================================================

sp=zeros(tf,3);         % setpoint trajectory

N=10;                   % prediction horizon
M=3;                    % control horizon

Q = diag([1 0.001 1]);  % state penalty
Pf = Q;                 % terminal state penalty
R = 0.01*eye(m);         % control penalty
%==========================================================================
% Simulation
%==========================================================================
    
% Simulation
xe_k_k1(:,1) = [zeros(n,1);zeros(nd,1)];
x07 = zeros(3,1);
x(:,1) = x07;

    for k = 1:tf
        
        %=============================
        % Calculate steady state target
        %=============================
        de = xe_k_k1(end-nd+1:end,k);
        xs_us = Mss*de;
        xs7 = xs_us(1:n);
        us7 = xs_us(end-m+1:end);
        
        %=============================
        % Solve the QP
        %=============================
        xe = xe_k_k1(1:n, k);
        dx = xe - xs7;
        [du,~] = CRHC2_57(A,B,N,M,Q,R,Pf,[],[],[],[],[],[],dx);
        u(:,k) = du(1:m) + us7;  
        
        %=============================
        % Update the observer state
        %=============================
        y(:,k) = C*x(:,k);
        xe_k_k = xe_k_k1(:,k) + Le*(y(:,k) - Ce*xe_k_k1(:,k));
        xe_k_k1(:,k+1) = Ae*xe_k_k + Be*u(:,k); %Updates here
        %=============================
        % Update the process state
        %=============================
        x(:,k+1) = A*x(:,k) + B*u(:,k) + Bp*d(k);
        
        %=============================        
        % Store current variables in log 
        %=============================

   end % simulation loop
 
        %%
%==========================================================================
% Plot results
%==========================================================================
        
   % plot the states, the state estimations, and the input and report them
   % in the report.
   if example == 'a'
   figure (1)
   subplot(1,3,1)
   plot(1:1:tf,x(1:n,1:end-1),'Linewidth',2);
   legend('c','T','h','Location','southeast');
   xlabel('Time'); ylabel('Measurement');
   title('$${x}$$','Interpreter','Latex');
   subplot(1,3,2)
   plot(1:1:tf,xe_k_k1(1:n,1:end-1),'Linewidth',2);
   legend('c_e','T_e','h_e','Location','southeast');
   xlabel('Time'); ylabel('Estimation');
   title('$$\hat{x}$$','Interpreter','Latex');
   subplot(1,3,3)
   plot(1:1:tf,u(1:m,1:end),'Linewidth',2);
   legend('u_1','u_2','Location','southeast');
   xlabel('Time'); ylabel('Inputs');
   title('$${u}$$','Interpreter','Latex');
   suptitle('System a');
   end
   if example == 'c'
   figure (1)
   subplot(1,3,1)
   plot(1:1:tf,x(1:n,1:end-1),'Linewidth',2);
   legend('c','T','h','Location','southeast');
   xlabel('Time'); ylabel('Measurement');
   title('$${x}$$','Interpreter','Latex');
   subplot(1,3,2)
   plot(1:1:tf,xe_k_k1(1:n,1:end-1),'Linewidth',2);
   legend('c_e','T_e','h_e','Location','southeast');
   xlabel('Time'); ylabel('Estimation');
   title('$$\hat{x}$$','Interpreter','Latex');
   subplot(1,3,3)
   plot(1:1:tf,u(1:m,1:end),'Linewidth',2);
   legend('u_1','u_2','Location','east');
   xlabel('Time'); ylabel('Inputs');
   title('$${u}$$','Interpreter','Latex');
   suptitle('System c');
   end