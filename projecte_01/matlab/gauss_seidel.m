function [x, it] = gauss_seidel(A, b, x0, itmax, tol)
%--------------------------------------------------------------------------
% Description: this function computes the solution to the linear system 
% A x = b using Gauss-Seidel method 
%--------------------------------------------------------------------------
% Inputs: 
%   - A         Linear system matrix    [n x n]
%   - b         Linear system vector    [n x 1]
%   - x0        Initial solution        [n x 1]
%   - itmax     Max number of iterations
%   - tol       Tolerance to check convergence
%--------------------------------------------------------------------------
% Outputs: 
%   - x         Solution to A x = b
%--------------------------------------------------------------------------

% Output variables
x = x0;         % Solution to the linear system
it = 0;         % Iteration counter

n = length(b);  % Linear system dimension
conv_cond = 0;  % Convergence condition boolean variable

% Gauss-Seidel method
while (it <= itmax && conv_cond == 0)
    % Increase iteration counter counter
    it = it + 1;
    % Compute next step
    x_temp = x; % Temporary vector, to check convergence condition later
    for k = 1:n
        x(k) = (b(k) - A(k,:)*x + A(k,k)*x(k))/A(k,k);
    end
    % Convergence condition, uses infinite norm
    if max(abs(x_temp-x)) < tol
        conv_cond = 1;
    end
end

end