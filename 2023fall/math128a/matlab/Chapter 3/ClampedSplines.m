function Splines = ClampedSplines(x,f,df)
%
% This code implements the clamped splines
%
% Written by Ming Gu for Math 128A, Fall 2008
% Updated by Ming Gu for Math 128A, Spring 2015
%
n = length(x);
h = diff(x(:));
rhs = 3 * diff([df(1);diff(f(:))./h;df(2)]);

A = diag(h,1)+diag(h,-1)+2*diag([[0;h]+[h;0]]);

%
% compute the c coefficients. This is a simple 
% but very slow way to do it.
%
c         = A \ rhs;
d         = (diff(c)./h)/3;
b         = diff(f(:))./h-(h/3).*(2*c(1:n-1)+c(2:n));
Splines.a = f(:);
Splines.b = b;
Splines.c = c;
Splines.d = d;
