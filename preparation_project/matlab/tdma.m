function [x] = tdma(A, b)
%--------------------------------------------------------------------------
% Description: this function computes the solution to the linear system 
% A x = b, where A is a tri-diagonal matrix 
%--------------------------------------------------------------------------
% Inputs: solution vectors
%   - A         [n x 3] matrix, where col 1 is lower-digonal, col 2 is 
%                the main diagonal and col 3 is the upper-diagonal
%   - b         Linear system vector    [n x 1]
%--------------------------------------------------------------------------
% Outputs: 
%   - x         Solution to A x = b
%--------------------------------------------------------------------------

n = length(b);  % Linear system dimension

% P coefficients vector
P = zeros(n, 1);
P(1) = A(1,3)/A(1,2);
for i = 2:n-1
    P(i) = A(i,3)/(A(i,2) - A(i,1)*P(i-1));
end

% R coefficients vector
R = zeros(n, 1);
R(1) = b(1)/A(1,2);
for i = 2:n
    R(i) = (b(i) + A(i,1)*R(i-1))/(A(i,2) - A(i,1)*P(i-1));
end

% Linear system solution
x = zeros(n, 1);
x(n) = R(n);
for i = n-1:-1:1
    x(i) = P(i)*x(i+1) + R(i);
end

end