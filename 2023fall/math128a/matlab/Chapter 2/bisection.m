function [x, out] = bisection(Fcn, Intv, params)
% 
%  To find a solution to f(x) = 0 given the continuous function
%  f on the interval [a,b], where f(a) and f(b) have
%  opposite signs:
% 
%  INPUT:  function f(x) defined by function handle Fcn, 
%          interval [a,b]= [Intv.a, Intv.b]
%          tolerence params.tol, max # of iterations = params.MaxIt
%  OUTPUT: root x, and data structure out. 
%          The success flag out.flg, is 0 for successful
%          execution and non-zero otherwise. out.it is the number
%          of iterations to reach within tolerance.
%
% Written by Ming Gu for Math128A, Spring 2022

TOL = params.tol;
NO  = params.MaxIt;
a    = Intv.a;
b    = Intv.b;
if (a > b)
   a = Intv.b;
   b = Intv.a;
end
fa   = sign(Fcn(a));
fb   = sign(Fcn(b));
if (fa*fb >0)
    error('Initial Interval may not contain root',msg);
end
if a==b
    error('Initial values for a and b must not equal',msg);
end

It = 0;
out.x =[a;b];
out.f =[Fcn(a);Fcn(b)];
while (It <= NO)
   c = (a+b)/2;
   out.x = [out.x;c];
   out.f =[out.f;Fcn(c)];
   fc = sign(Fcn(c));
   if (fc ==0)
      x = c;
      out.flg = 0;
      out.it  = It;
      return;
   end
   if (fc * fa < 0)
      b  = c;
      fb = fc;
   else
      a  = c;
      fa = fc;
   end
   if (abs(b-a)<=TOL)
      x = (a+b)/2;
      out.flg = 0;
      out.it  = It;
      return;
   end
   It = It + 1;
end
out.flg =1;
out.it = NO;
x = (a+b)/2;
