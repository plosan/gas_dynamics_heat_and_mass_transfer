clear;
close all;
clc;
% 
% v = [1,3,5];
% A = rand(5,5)
% 
% A(v,:) = []
% A(:,v) = []


v = randi(10,8,1)

v(end-4:end)

v([2:2:8]) = 15;

v



% A = [
%     3   1   -1;
%     1   4  2;
%     3   -2  6;
% ];
% b = rand(3,1);
% cond(A)
% 
% x_matlab = linsolve(A,b)
% 
% x0 = x_matlab + 0.5 ;
% itmax = 50;
% tol = 1e-16;
% 
% [x_gs, it] = gaussSeidel(A, b, x0, itmax, tol);
% 
% res = A*x_gs - b;
% norm(res)
% it