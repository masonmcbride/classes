%
% this code demonstrates how the different ODE methods work. 
%
% Written by Ming Gu for Math 128a, Fall 2021
% Updated Fall 2022
%
format short e;
Intv     = [1;3/2];
tol      = 1e-6;
stepsize = [1e-4;1e-1];
alpha    = -1;
N        = 200;
%
%  the explicit 4th order Runge-Kutta method 
%
figure(1);
[w_RK4,t_RK4] = RK4(@ty3, Intv, alpha, N);
plot(t_RK4,w_RK4);
hold on;
plot(t_RK4,w_RK4,'rd');
ff1 = ['Runge Kutta 4th order'];
ff2 = ['Number of mesh points = ', num2str(length(t_RK4))];
ff  = strvcat(ff1,ff2);
title(ff,'FontSize',18);
%
% the Runge-Kutta-Fehlberg method 
%
figure(2);
[w_RKF,t_RKF] = RKF(@ty3,Intv,alpha,tol,stepsize);
plot(t_RKF,w_RKF);
hold on;
plot(t_RKF,w_RKF,'rd');
ff1 = ['Runge-Kutta-Fehlberg method'];
ff2 = ['Number of mesh points = ', num2str(length(t_RKF))];
ff  = strvcat(ff1,ff2);
title(ff,'FontSize',18);
%
% Adams 4th order Predictor-Corrector 
%
figure(3);
[w_adams4,t_adams4] = Adams4PC(@ty3, Intv, alpha, N);
plot(t_adams4,w_adams4);
hold on;
plot(t_adams4,w_adams4,'rd');
ff1 = ['Adams 4th order PC'];
ff2 = ['Number of mesh points = ', num2str(length(t_adams4))];
ff  = strvcat(ff1,ff2);
title(ff,'FontSize',18);
%
% Adams Variable Step-size PC method
%
figure(4);
[w_adamsadpt,t_adamsadpt] = AdamsAdaptPC(@ty3,Intv,alpha,tol,stepsize);
plot(t_adamsadpt,w_adamsadpt);
hold on;
plot(t_adamsadpt,w_adamsadpt,'rd');
ff1 = ['Adams Variable Step-size PC method'];
ff2 = ['Number of mesh points = ', num2str(length(t_adamsadpt))];
ff  = strvcat(ff1,ff2);
title(ff,'FontSize',18);

