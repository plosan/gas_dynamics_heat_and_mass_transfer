function [T, it, res] = modifiedGaussSeidel(A, b, T0, N, L, itmax, tol)
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
T = T0;         % Solution to the linear system
it = 0;         % Iteration counter
res = 1;        % Residue

% n = (N+2)*(L+2);  % Linear system dimension
conv_cond = 0;      % Convergence condition boolean variable

while(it <= itmax && conv_cond == 0)    
    % Increase iteration counter
    it = it + 1;
    % Auxiliary vector
    T_aux = T; % Temporary vector, to check convergence condition later
    % Compute next iteration
    
    % ROW j=1 (LOWER WALL NODES)
    for i = 2:N+1
        T(i) = (b(i) + A(i,2)*T(i-1) + A(i,4)*T(i+1) + A(i,5)*T(i+(N+2)))/A(i,3);
    end
    % ROW 2<=j<=L+1 (MID ROWS)
    for j = 2:L+1
        % LEFT WALL NODE
        id = 1 + (j-1)*(N+2);
        T(id) = (b(id) + A(id,1)*T(id-(N+2)) + A(id,4)*T(id+1) + A(id,5)*T(id+(N+2)))/A(id,3);
        % INTERNAL NODES
        for i = 2:N+1
            id = i + (j-1)*(N+2);
            T(id) = (b(id) + A(id,1)*T(id-(N+2)) + A(id,2)*T(id-1) + A(id,4)*T(id+1) + A(id,5)*T(id+(N+2)))/A(id,3);
        end
        % RIGHT WALL NODE
        id = j*(N+2);
        T(id) = (b(id) + A(id,1)*T(id-(N+2)) + A(id,2)*T(id-1) + A(id,5)*T(id+(N+2)))/A(id,3);
    end
    
    % ROW j=L+2 (UPPER WALL NODES)
    for i = 2:N+1
        id = i + (L+1)*(N+2);
        T(id) = (b(id) + A(id,1)*T(id-(N+2)) + A(id,2)*T(id-1) + A(id,4)*T(id+1))/A(id,3);
    end
    
    % Convergence condition, uses infinite norm
    res = max(abs(T_aux-T));
    if res < tol
        conv_cond = 1;
    end
end


end