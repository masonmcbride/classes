function F = NDD1(x,f)
%
% This function implements Newton's Divided Difference Algorithm
% F is the vector of coefficients
%
% Updated by Ming Gu for Math 128A, Spring 2021
%
N = length(x);
F = f;
for k=2:N
    for j = N:-1:k
       F(j) = (F(j)-F(j-1))/(x(j)-x(j-k+1));
    end
end


   
