function [Int,flg, fcnt,level] = AdaptSimpson(Fun,interv,tol,L)
% 
% Driver routine to Adaptive Simpson's Rule code
%
% On input: 
%   Fun is the name of function to be integrated
%   interv is the interval to be integrated over
%   tol is tolerance
%   L is the maximum level of Adaptive Simpson recursion
%
% On output
%   Int is the integral
%   flg success flag: flg = 0 indicates success; 
%      flg == 1 indicates maximum level recursion reached.
%   fcnt is total number of function evaluations.
%   level is total number of recursions performed
%
% This routine calls the Adaptive Simpson recursion code 
% AdaptSimpson2 to perform the actual integration.
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
if (tol < eps | L < 1)
    error('Invalid Tolerance or Level');
end

%
% Evaluate the function at three nodal points
%
x = [a;(a + b)/2;b];
try
   f = Fun(x);
catch
   error('InvalidFunctionSupplied');
end
fx      = [x, f];
simpson = ([1 4 1] * f)*(b-a)/6;
[Int,flg,fcnt,level] = AdaptSimpson2(Fun,tol,L,fx,simpson);
fcnt    = fcnt + 3;
level   = L - level + 1;

function [Int,flg,fcnt,level] = AdaptSimpson2(Fun,tol,L,fx,simpson)
% 
% Recursive Adaptive Simpson's Rule
% 
% Written by Ming Gu for Math 128A, Fall 2022
% 
% Evaluate the function at three nodal points
%
xnew = [fx(1,1)+fx(2,1);fx(2,1)+fx(3,1)]/2;
fnew = Fun(xnew);
fcnt = 2;
h2   = (fx(3,1)-fx(1,1))/4;
simpson1 = sum([1 4 1] .* [fx(1,2) fnew(1) fx(2,2)])*h2/3;
simpson2 = sum([1 4 1] .* [fx(2,2) fnew(2) fx(3,2)])*h2/3;
Int      = simpson1+simpson2;
level    = L;
if (abs(Int-simpson)<15*tol)
    flg      = 0;
    return;
end
if (L == 1)
    flg = 1;
    return;
end
fx1  = [fx(1,1) fx(1,2);xnew(1) fnew(1); fx(2,1) fx(2,2)];
fx2  = [fx(2,1) fx(2,2);xnew(2) fnew(2); fx(3,1) fx(3,2)];
[Int1,flg1,fcnt1,level1] = AdaptSimpson2(Fun,tol/2,L-1,fx1,simpson1);
[Int2,flg2,fcnt2,level2] = AdaptSimpson2(Fun,tol/2,L-1,fx2,simpson2);
Int  = Int1 + Int2;
flg  = max(flg1, flg2);
fcnt = 2 + fcnt1 + fcnt2;
level= min(level1,level2);

