function Int = Simpson(Fun,interv)
% 
% Simpson's Rule
% 
% On input: 
%   Fun is the name of function to be integrated
%   interv is the interval to be integrated over
%
% On output
%   Int is the integral
%
% Written by Ming Gu for Math 128A, Fall 2022
% 
% Validate input arguments
%
if (length(interv)~=2)
    error('Invalid Interval');
end
a = interv(1);
b = interv(2);
if (a>b)
    error('Invalid Interval');
end
if (a == b)
    Int = 0;
    return;
end
%
% Evaluate the function at three nodal points
%
h = (b - a)/2;
c = (a + b)/2;
x = [a;c;b];
try
   fx = Fun(x);
catch
   error('InvalidFunctionSupplied');
end
Int = h*(fx(1)+fx(3)+4*fx(2))/3;







