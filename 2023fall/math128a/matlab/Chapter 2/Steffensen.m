function [x, out] = SteffensenMethod(Fcn, x0, params)
% 
%  To find a solution to x = f(x) given the continuous function
%  f and initial guess x0
%  opposite signs:
% 
%  INPUT:  function f(x) defined by function handle Fcn, 
%          initial guess x0, 
%           tolerence params.tol, max # of iterations = params.MaxIt
%  OUTPUT: root x, and data structure out. 
%          The success flag out.flg, is 0 for successful
%          execution and non-zero otherwise. out.it is the number
%          of iterations to reach within tolerance.
%
% Written by Ming Gu for Math128A, Spring 2021

TOL = params.tol;
NO  = params.MaxIt;

It = 1;
x  = x0;
out.x = [x];
out.f = [Fcn(x)];
while (It <= NO)
   y = Fcn(x);
   z = Fcn(y);
   p = x - ((y-x)^2)/((z-y)-(y-x));
   out.x = [out.x;x];
   out.f =[out.f;y;z];
   if (abs(p-x) <=TOL)
      x = p;
      out.flg = 0;
      out.it  = It;
      out.x = [out.x;x];
      return;
   end
   x  = p;
   It = It + 1;
end
out.flg =1;
out.it = NO;
