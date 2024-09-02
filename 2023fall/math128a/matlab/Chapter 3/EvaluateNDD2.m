function f = EvaluateNDD2(xnew,x,F)
%
% This function evaluates the Hermite interpolating polynomial given
% point xnew, nodes x, and NDF coefficients F
%
n = length(x);
m = length(xnew);
f = F(2*n)*ones(m,1);
z = kron(x(:),ones(2,1));
for k=2*n-1:-1:1
   f = F(k) + f .* (xnew-z(k));
end

