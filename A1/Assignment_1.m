clc
clear all

%Q1
Ac = [0 1; 0.5 0]; Bc = [0;1];
Cc = [1 0]; Dc = 0;
sysc = ss(Ac,Bc,Cc,Dc); %% System (2)
h = 0.1;
sysd = c2d(sysc,h); %% System (3)
A3 = expm(Ac*h);
fun0 = @(t) expm(Ac*t)*Bc;
B3 = integral(fun0,0,h,'ArrayValued',true);
syms k;
fun = @(t) expm(Ac*t)*Bc;
m = (integral(fun,0,0.05,'ArrayValued',true));
n = (integral(fun,0,0.05,'ArrayValued',true));
Aa = [expm(Ac*h) m; 0 0 0;];
Ba = [n; 1]; Ca = [1 0 0];
sysd2=ss(Aa,Ba,Ca,Dc); %%System (4)

%Q3
obs1 = obsv(sysc); rank1 = rank(obs1);
obs2 = obsv(sysd); rank2 = rank(obs2);
obs3 = obsv(sysd2); rank3 = rank(obs3);
%C_t = [1 0];
%C_t1 = [2 1];
%C_t1*Ac

x2 = 1; x1 = sqrt(0.5)*x2;
C_t1 = [x1 x2];
obs_t = [C_t1; C_t1*Ac];
r_obs_t = rank(obs_t)
C_t12 = [x1 x2 0];
obs_td = [C_t1; C_t1*A3];
r_obs_td = rank(obs_td);

%Q5
sysc5 = ss(Ac,Bc,C_t1,Dc);
sysd5 = c2d(sysc5,h);

%Q6
p_c = [-4+6*1i, -4-6*1i];
Kc = place(Ac,Bc,p_c);
Ad6 = expm((Ac-Bc*Kc)*h);
p_d = eig(Ad6);
Kd = place(A3,B3,p_d);
sys3_step = ss(A3-B3*Kd, B3, Cc,0);
Tfinal = 20;

Kd4_0 = [Kd 0];
sys4_step = ss(Aa-Ba*Kd4_0, Ba, Ca, 0);
figure (1)
subplot(2,1,1)
step(sys3_step,Tfinal)
title('System 3');
subplot(2,1,2)
step(sys4_step, Tfinal)
title ('System 4');

Kd4_delay = place(Aa,Ba,[p_d ;0]);
sys4del = ss((Aa-Ba*Kd4_delay),Ba,Ca,0);

figure(2)
subplot(2,1,1)
step(sys4_step, Tfinal)
title('System 4 - Old Feedback');
subplot(2,1,2)
step(sys4del,Tfinal)
title ('System 4 - New Feedback');







