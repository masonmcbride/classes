function [y,out] = Neville(xx,x,Q)
%
% Neville's method
% 
% inputs:
%    n = order of interpolation (n+1 = # of points)
%    x(1),...,x(n+1)    x coords
%    Q(1),...,Q(n+1)    function values at x
%    xx=evaluation point for interpolating polynomial p
%
%
% On output
%    y     = p(xx): 
%    out.Q = intermediate Q values.
%    out.Q(j): approximation by interpolating at points x(j) ... x(n+1)
%
% Written by Ming Gu for Math 128A, Spring 2021 
% 

N = length(x);
n = N - 1;

for i = n:-1:1
   for j = 1:i
      Q(j) = (xx-x(j))*Q(j+1) - (xx-x(j+N-i))*Q(j);
      Q(j) = Q(j)/(x(j+N-i)-x(j));
   end
end

y     = Q(1);
out.Q = Q;

