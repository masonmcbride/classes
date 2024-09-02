function Int = CompSimpson(Fun,n,interval)
% 
% Composite Simpson's Rule
% 
% On input: 
%   FunIn is the name of function to be integrated
%   interv is the interval to be integrated over
%   n is the number of nodal points to be used.
%
% On output
%   Int is the integral
%
% Written by Ming Gu for Math 128A, Fall 2022
% 
% Validate input arguments
%
if (n < 3 | mod(n,2)==0)
    error('Invalid Number of Nodes');
end
if (length(interval)~=2)
    error('Invalid Interval');
end
a = interval(1);
b = interval(2);
if (a>b)
    error('Invalid Interval');
end
if (a == b)
    Int = 0;
    return;
end
%
% Evaluate the function at nodal points
%
h = (b - a)/(n-1);
x = a + (0:n-1)'*h;
try
   fx = Fun(x);
catch
   error('InvalidFunctionSupplied');
end
Int = h*(fx(1)+fx(n)+4*sum(fx(2*(1:(n-1)/2)))+2*sum(fx(2*(2:(n-1)/2)-1)))/3;

