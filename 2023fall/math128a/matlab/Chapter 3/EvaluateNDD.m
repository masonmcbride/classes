function f = EvaluateNDD(xnew,x,F)
%
% This function evaluates the interpolating polynomial given
% point xnew, nodes x, and NDF coefficients F
%
n = length(x);
m = length(xnew);
f = F(n)*ones(m,1);
for k=n-1:-1:1
   f = F(k) + f .* (xnew-x(k));
end

