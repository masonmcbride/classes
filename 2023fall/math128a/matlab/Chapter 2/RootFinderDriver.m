%
% This is a driver for three root finders

% Written by Ming Gu for Math 128A
%

format long g
params.tol   = 1e-15;
params.MaxIt = 100;
close all;
%
% bisection method
%
Intv.a = 2.5;
%%Intv.b = 10;
Intv.b = 3.5;
[x_bisect, out_bisect]= bisection(@myfunc,Intv,params);

%
% fixed point
%
x0= Intv.a;
[x_fixedpoint, out_fixedpoint] = fixedpoint(@myfunc1, x0, params);

%
% Secant method
%
x0 = 2.5;
x1 = 4;
[x_secant, out_secant]= SecantMethod(@myfunc,x0,x1,params);

%
% Newton's Method
%
x0= Intv.b;
x0= Intv.a;
%%x0= 2.5; %% method will quickly converge with this choice of x0 
[func,dfunc,x_Newton, out_Newton]=NewtonMethod(@myfunc, @dmyfunc, x0, params);
%
% Report results
%
fprintf('\n');
if (out_bisect.flg==0)
   disp(['Bisection succeeded.']);
   disp(['# of iterations = ',num2str(out_bisect.it), '; root = ', num2str(x_bisect(end))]);
   fx = myfunc(x_bisect(end));
   disp(['Residual Error = ', num2str(fx)]);
   fprintf('\n');
   fprintf('\n');

   figure(1);
   grid on;
   xx = min(out_bisect.x)+((0:1000)/1000)*(max(out_bisect.x)-min(out_bisect.x));
   fxx= myfunc(xx);
   plot(xx,fxx,'b.-');
   hold on
   [xx,junk] = sort(out_bisect.x);
   plot(xx,out_bisect.f(junk),'kd');
   title(['Bisection Method with ', num2str(length(xx)), ' Points']);
 
else
   disp(['Bisection failed.']);
   fprintf('\n');
   fprintf('\n');
end
if (out_fixedpoint.flg==0)
   disp(['Fixedpoint succeeded.']);
   disp(['# of iterations = ',num2str(out_fixedpoint.it), '; root = ', num2str(x_fixedpoint(end))]);
   fx = myfunc(x_fixedpoint(end));
   disp(['Residual Error = ', num2str(fx)]);
   fprintf('\n');
   fprintf('\n');

   figure(2);
   grid on;
   xx = min(out_fixedpoint.x)+((0:1000)/1000)*(max(out_fixedpoint.x)-min(out_fixedpoint.x));
   fxx= myfunc(xx);
   plot(xx,fxx,'b.-');
   hold on
   [xx,junk] = sort(out_fixedpoint.x);
   plot(xx,myfunc(xx),'kd');
   title(['FixedPoint Method with ', num2str(length(xx)), ' Points']);
else
   disp(['Fixedpoint failed.']);
   fprintf('\n');
   fprintf('\n');
end

if (out_secant.flg==0)
   disp(['Secant Method succeeded.']);
   disp(['# of iterations = ',num2str(out_secant.it), '; root = ', num2str(x_secant(end))]);
   fx = myfunc(x_secant(end));
   disp(['Residual Error = ', num2str(fx)]);
   fprintf('\n');
   fprintf('\n');

   figure(3);
   grid on;
   xx = min(out_secant.x)+((0:1000)/1000)*(max(out_secant.x)-min(out_secant.x));
   fxx= myfunc(xx);
   plot(xx,fxx,'b.-');
   hold on
   [xx,junk] = sort(out_secant.x);
   plot(xx,out_secant.f(junk),'kd');
   title(['SecantMethod with ', num2str(length(xx)), ' Points']);
else
   disp(['Secant failed.']);
end

if (out_Newton.flg==0)
   disp(['Newton succeeded.']);
   disp(['# of iterations = ',num2str(out_Newton.it), '; root = ', num2str(x_Newton(end))]);
   fx = myfunc(x_Newton(end));
   disp(['Residual Error = ', num2str(fx)]);
   fprintf('\n');
   fprintf('\n');

   figure(4);
   grid on;
   xx = min(out_Newton.x)+((0:1000)/1000)*(max(out_Newton.x)-min(out_Newton.x));
   fxx= myfunc(xx);
   plot(xx,fxx,'b.-');
   hold on
   [xx,junk] = sort(out_Newton.x);
   plot(xx,out_Newton.f(junk),'kd');
   title(['Newton Method with ', num2str(length(xx)), ' Points']);

else
   disp(['Newton failed.']);
end

