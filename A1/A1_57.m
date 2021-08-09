clear;clc;close all;
%%You are supposed to fill this template. Remeber not to change variable
%%names! When you run this file, you should see all these variables in the Workspace with
%%the right values! 
%% Discrete Model
Ac=[0 1;0.5 0];Bc=[0;1];Cc=[1 0];Dc=0;h=0.1;
sysc = ss(Ac,Bc,Cc,Dc); %System 2 - Continuous
%Use given parameters and find A, B, and C in Question #1

A= expm(Ac*h);
fun0 = @(t) expm(Ac*t)*Bc;
B= integral(fun0,0,h,'ArrayValued',true);
C= Cc; D = Dc;

sysd = c2d(sysc,h); %System 3 - Discrete

%% Delayed model
% define Aa, Ba, Ca, Da according to Question #2
fun = @(t) expm(Ac*t)*Bc;
B1 = expm(Ac*0.05)*(integral(fun,0,0.05,'ArrayValued',true));
B2 = (integral(fun,0,0.05,'ArrayValued',true));

Aa= [expm(Ac*h) B1; 0 0 0];
Ba= [B2; 1];
Ca= [Cc 0];
Da= 0;
sysd2 = ss(Aa,Ba,Ca,Da,h); %System 4 - Delayed system
%% Controllability and Observability
%find the rank of controllablity and observability matrices according to
%Question #3


% the following parameters should be the ranks of the observability and 
% the controllability matrices, accordingly


ctr_sys2 = ctrb(sysc); sys2_c= rank(ctr_sys2);
obs_sys2 = obsv(sysc); sys2_o= rank(obs_sys2); 
ctr_sys3 = ctrb(sysd); sys3_c= rank(ctr_sys3);
obs_sys3 = obsv(sysd); sys3_o= rank(obs_sys3);
ctr_sys4 = ctrb(sysd2); sys4_c= rank(ctr_sys4);
obs_sys4 = obsv(sysd2); sys4_o= rank(obs_sys4);

%Q4 

x2 = 1; x1 = sqrt(0.5)*x2; C_unobs = [x1 x2];
unobs_2 = [C_unobs; C_unobs*Ac]; unobs_3 = [C_unobs; C_unobs*A]; 
rank_sys2_obs = rank(unobs_2); rank_sys3_obs = rank(unobs_3); 

%Q6
%% Controller design
landa1=-4+6*1i;
landa2=-4-6*1i;
Kc = place(Ac,Bc,[landa1, landa2]);
Ac_k = (Ac - Bc*Kc);
Ad_k = expm(Ac_k*h);
fun_k = @(t) expm(Ac_k*t)*Bc;
Bd_k = integral(fun_k,0,h,'ArrayValued',true);
pol_sys3_K = eig(Ad_k);
%calculate desired poles for the discrete time system (3) and define them as p1 and
%p2 as it is asked in Question #6
p1= pol_sys3_K(1,1);
p2= pol_sys3_K(2,1);

%define the feedback gain for the discrete time system (3) as K1
K1= place(A,B,[p1,p2]);
K1_sys4 = [K1 0];
sys3K = ss(A-B*K1,B,C,0,h);
sys4K = ss((Aa-Ba*K1_sys4),Ba,Ca,Da,h);
Tf = 10;
figure (1)
step(sys3K,0:h:Tf)
title ('Step response of system 3 with controller K3');

y1 = step(sys3K,0:h:Tf);
figure (2)
step(sys4K,0:h:Tf)
title ('Step response of system 4 with controller K3');
y2 = step(sys4K,0:h:Tf);

%define the feedback gain for the delayed discrete time system (4) as K2
K2= place(Aa,Ba,[p1,p2, 0]);
sys4K2 = ss((Aa-Ba*K2),Ba,Ca,Da,h);  
figure (3)
step(sys4K2,0:h:Tf)
y3 = step(sys4K2,0:h:Tf);

%plot the step response of the systems in one figure. Your figure should
%have labels and legend.

t = 0:h:10;
figure (4)
plot(t,y1,t,y2,t,y3);
title ('Response of system 3 and 7 with previous and modified controller K')
legend ('Sys 3 - K3','Sys 7 - K3','Sys 7 - K7'); 
xlabel ('Time (s)');
ylabel ('y(k) (radians)');

%% Steady State
ys=pi/6;
AdK2 = Aa-Ba*K2;

x1s = ys; syms x2s x3s us;
eqn1 = x1s == AdK2(1,1)*x1s + AdK2(1,2)*x2s + AdK2(1,3)*x3s + Ba(1)*us;
eqn2 = x2s == AdK2(2,1)*x1s + AdK2(2,2)*x2s + AdK2(2,3)*x3s + Ba(2)*us;
eqn3 = x3s == AdK2(3,1)*x1s + AdK2(3,2)*x2s + AdK2(3,3)*x3s + Ba(3)*us;
eqn = [eqn1; eqn2; eqn3];
[x2ss, x3ss, uss] = solve(eqn,x2s,x3s,us);
steady_u = uss;
steady_x = [x1s; x2ss; x3ss];

% plot the system output as explained in Question #7. Your figure should
%have labels and legend.
t = 1:h:5;
x = zeros(3,length(t));

for k = 1:length(t)
   x(:,k+1) = (Aa-Ba*K2)*x(:,k) + Ba*steady_u;
    y(:,k) = Ca*x(:,k);
end
figure (5)
 plot(t,y);
 title('Steady State Response for input Us = 17.94');
 xlabel ('Time (s)');
 ylabel ('y(k) (radians)');



    


%% disturbance
Bd=[0;1;0];
% define Ae, Be, Ce, and De as asked in Question #8
Ae= [Aa,Bd; 0,0,0,1];
Be= [B2; 1; 0];
Ce= [1,0,0,0];
De= 0;

sys6_poles = eig(Ae);
%define the rank of the controllability and observability matrices as

sys6_c= [Be Ae*Be Ae*Ae*Be Ae*Ae*Ae*Be];
sys6_o= [Ce; Ce*Ae; Ce*Ae*Ae; Ce*Ae*Ae*Ae];
r_6c = rank(sys6_c);
r_60 = rank(sys6_o); % Detectable but not stable.

%% Feedback gain
%define the controller gain as K3 according to Question #9
K3= place(Ae, Be, [p1,p2,0,1]);




%% Observer
%define L as the observer gain according to Question #10
L = place(Ae',Ce',[0.1,0.2,0.3,0.4]);

