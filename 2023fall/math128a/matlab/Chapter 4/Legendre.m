function [c,x] = Legendre(n)
%
% this routine computes the weights and nodes 
% for the n-node Gaussian quadrature 
% through the n-th order Legendre polynomial.
%
% Written by Ming Gu for Math 128A, Fall 2008
%
b = (1:n-1)';
b = b./sqrt((2*b-1).*(2*b+1));
B = diag(b,1)+diag(b,-1);
[Q,D] = eig(B);
x = diag(D);
x(abs(x)<1e-15) = 0;
c = 2*(Q(1,:).^2)';
