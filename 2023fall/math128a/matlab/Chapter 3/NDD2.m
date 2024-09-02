function F = NDD2(x,f,df)
%
% This function implements Newton's Divided Difference Algorithm
% for Hermite Interpolation. f is the vector of function values
% and df vector of derivatives.
%
% Ming Gu for Math 128A, Spring 2021
%
N = length(x);
x = x(:);
xx = reshape(repmat(x',2,1),2*N,1);
f = f(:);
df= df(:);
F  = reshape(repmat(f',2,1),2*N,1);
NN = N * 2;
F(2*(1:N))     = df;
F(1+2*(1:N-1)) = (f(2:N)-f(1:N-1))./(x(2:N)-x(1:N-1));
for k=3:2*N
    for j = 2*N:-1:k
       F(j) = (F(j)-F(j-1))/(xx(j)-xx(j-k+1));
    end
end


   
