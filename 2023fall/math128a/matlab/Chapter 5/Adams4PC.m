function [w,t] = Adams4PC(Fcn, Intv, alpha, N)
% On input: 
%   Fcn is the name of function to be integrated
%   interv is the interval to be integrated over
% 
%   The problem: y' = f(t,y), y(a) = alpha, a<= t <= b
%   where Intv = [a b], and a call to Fcn with 
%   argument (t, y) returns f(t,y).
%
%   Adams4PC uses the Adams 4th order Predictor-Corrector method 
%   to approximate y at N+1 equi-spaced points on [a, b].
%
% On output
%   t contains the (equi-spaced) mesh points, w the
%   function values at these points. 
%
% Written by Ming Gu for Math 128A, Fall 2008
% Updated Spring 2021
% 
a    = Intv(1);
b    = Intv(2);
h    = (b-a)/N;
w    = zeros(N+1,1);
t    = a + h*(0:N)';
w(1) = alpha;
%
% RK4 for the first 3 steps
%
h2   = h/2;
for i = 1:3
    k1 = h* Fcn(t(i),w(i));
    k2 = h* Fcn(t(i)+h2,w(i)+k1/2);
    k3 = h* Fcn(t(i)+h2,w(i)+k2/2);
    k4 = h* Fcn(t(i)+h,w(i)+k3);
    w(i+1) = w(i) + (k1+2*k2+2*k3+k4)/6;
end
%
% main loop
%
p = h*[-9/24  37/24 -59/24 55/24];
c = h*[ 1/24 -5/24   19/24 9/24 ];
f = Fcn(t(1:4), w(1:4)); 
for i = 4:N
    wp     = w(i) + p*f;
    fp     = Fcn(t(i+1),wp);
    w(i+1) = w(i) + c *[f(2:end);fp];
    f      =[f(2:end); Fcn(t(i+1),w(i+1))];
end

