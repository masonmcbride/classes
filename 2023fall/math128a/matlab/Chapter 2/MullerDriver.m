%
% this code finds a root of an input polynomial with
% Muller's method
%
% Written by Ming Gu for Math 128A, Fall 2022
%
n    = 3;
x_in = randn(n,1);
params.tolx   = 1e-8;
params.tolf   = 1e-8;
params.MaxIt = 100;
[fun, x, out] = Muller1(@mullerfun, x_in, params);
