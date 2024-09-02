function [x, out] = SecantMethod(Func, x0, x1, params)
% On input: 
%   Fun is the name of function 
%   x0, x1 are initial guesses, and tol tolerance.
%
% On output
%   fun contains all function values 
%   computed by Secant
%   out.flg = 0 means success; otherwise method failed.
%   x(end) is the root if out.flg = 0.
%   
% Written by Ming Gu for Math 128B, Spring 2010
% Updated by Ming Gu for Math 128A, Fall 2022
% 

N = params.MaxIt;
tol = params.tol;

out.flg = 1;
x(1)= x0;
x(2)= x1;
fun(1) = Func(x(1));
fun(2) = Func(x(2));
out.x = [x(1);x(2)];
out.f = [fun(1);fun(2)];

for k = 2:N
    if (fun(k) == fun(k-1))
       x = x(end);
       out.it = k;
       return;
    end
    x(k+1)   = x(k) - fun(k)/(fun(k)-fun(k-1))*(x(k)-x(k-1));
    fun(k+1) = Func(x(k+1));
    out.x    = [out.x;x(k+1)];
    out.f    = [out.f;fun(k+1)];
    if (abs(fun(k+1)) < tol)
       out.flg = 0;
       x = x(end);
       out.it = k;
       return;
    end
end


       

