function [Delta2, out] = Aitken(x)
% On input: 
%   linearly convergent finite sequence x
%
% On output
%   Aitken's Delta^2 Method.
%   out.fig = 0: success 
%   out.fig = 1: failure
%
% Written by Ming Gu for Math 128A, Spring 2015
% 
out.flg = 0;
delx    = x(2:end)-x(1:end-1);
del2x   = delx(2:end)-delx(1:end-1);
Delta2  = x(1:end-2)-(delx(1:end-1).^2)./del2x;
if(~isfinite(sum(Delta2)))
  out.flg= 1;
end

