function [A, b, d] = computeConstantDiscretizationCoefficients2(N, L, ...
    t_step, beta, alpha_g, Tg, dx, dy, Sx, Sy, V, rho, cp, lambda_w, lambda_e, lambda_s, lambda_n)
%--------------------------------------------------------------------------
% Description: this function computes the constant discretization
% coefficients of the problem.
%--------------------------------------------------------------------------
% Inputs: 
%   - N         Number of CVs in x direction
%   - L         Number of CVs in y direction
%   - t_step    Time step [s]
%   - beta      Time integration constant
%   - alpha_g   Left wall fluid heat transfer coefficient
%   - Tg        Left wall fluid temperature
%   - dx        Distances between nodes in x direction. dx(i) is the distance from node i to i+1
%   - dy        Distances between nodes in y direction. dy(i) is the distance from node i to i+1
%   - Sx        Surfaces in x axis
%   - Sy        Surfaces in y axis
%   - V         Nodes volume
%   - rho       (L1+L2+2)x(N1+N2+2) matrix containing the density of each node's material
%   - cp        (L1+L2+2)x(N1+N2+2) matrix containing the specific heat of each node's material
%   - lambda_w  Thermal conductivities matrix in west face
%   - lambda_e  Thermal conductivities matrix in east face
%   - lambda_s  Thermal conductivities matrix in south face
%   - lambda_n  Thermal conductivities matrix in northface
%--------------------------------------------------------------------------
% Outputs: 
%   - A         Matrix containing coefficients aS (column 1), aW (column
%               2), aP (column 3), aE (column 4), aN (column 5). Each row
%               corresponds to one node (see report). Size ((N+2)*(L+2)) x
%               ((N+2)*(L+2)).
%   - b         Vector containing independent terms. Size ((N+2)*(L+2)) x 1
%   - d         Vector containing rho*V*cp/t_step for each internal node.
%               Created for convenience.
%--------------------------------------------------------------------------


% MATRIX OF DISCRETIZATION COEFFICIENTS
A = zeros((L+2)*(N+2),5);   % 
d = zeros((L+2)*(N+2),1);
b = zeros((L+2)*(N+2),1);

for j = 2:L+1
    for i = 2:N+1
        % First version of the discretization coefficients
        aW = beta*lambda_w(j,i)*Sy(j)/dx(i-1);
        aE = beta*lambda_e(j,i)*Sy(j)/dx(i);
        aS = beta*lambda_s(j,i)*Sx(i)/dy(j-1);
        aN = beta*lambda_n(j,i)*Sx(i)/dy(j);
        dP = rho(j,i)*V(j,i)*cp(j,i)/t_step;
        aP = dP+aW+aE+aS+aN;
        % Node P ID
        id = i+(j-1)*(N+2);
        % Insert coefficients
        A(id,1) = aS;
        A(id,2) = aW;
        A(id,3) = aP;
        A(id,4) = aE;
        A(id,5) = aN;
        d(id,1) = dP; 
    end
end

% LOWER AND UPPER WALL NODES
for i = 2:N+1
    % LOWER WALL NODES
    id = i;
    A(id,5) = lambda_n(1,i)/dx(1);
    A(id,3) = A(id,5) + alpha_g;
    b(id) = alpha_g*Tg;
    % UPPER WALL NODES
    id = i+(L+1)*(N+2);
    A(id,1) = lambda_s(L+2,i)/dy(L+1);
    A(id,3) = A(id,1) + alpha_g;
    b(id) = alpha_g*Tg;
end

% LEFT AND RIGHT WALL NODES
for j = 2:L+1
    % LEFT WALL NODES
    id = 1+(j-1)*(N+2);
    A(id,4) = lambda_e(j,1)/dx(1);
    A(id,3) = A(id,4) + alpha_g;
    b(id) = alpha_g*Tg;
    % RIGHT WALL NODES
    id = j*(N+2);
    A(id,2) = lambda_w(j,N+2)/dx(N+1);
    A(id,3) = A(id,2) + alpha_g;
    b(id) = alpha_g*Tg;    
end

end