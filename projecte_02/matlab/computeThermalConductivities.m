function [lambda_w, lambda_e, lambda_s, lambda_n] = ...
    computeThermalConductivities(N1, N2, L1, L2, L3, x, y, x_cv, y_cv, dx, dy, lambda)
%--------------------------------------------------------------------------
% Description: this function computes the thermal conductivities on the
% control volume faces
%--------------------------------------------------------------------------
% Inputs: 

%   - N1        Number of CV in x direction for M1 and M3
%   - N2        Number of CV in x direction for M2 and M4
%   - L1        Number of CV in y direction for M1 and part of M2
%   - L2        Number of CV in y direction for part of M2 and part of M3
%   - L3        Number of CV in y direction for part of M3 and M4
%   - x         Nodes position in x direction
%   - y         Nodes position in y direction
%   - x_cv      CV walls position in x direction
%   - y_cv      CV walls position in y direction
%   - dx        Distances between nodes in x direction. dx(i) is the distance from node i to i+1
%   - dy        Distances between nodes in y direction. dy(i) is the distance from node i to i+1
%   - lambda    (L1+L2+2)x(N1+N2+2) matrix containing the thermal conductivity of each node's material
%--------------------------------------------------------------------------
% Outputs: 
%   - lambda_w  Thermal conductivities matrix in west face
%   - lambda_e  Thermal conductivities matrix in east face
%   - lambda_s  Thermal conductivities matrix in south face
%   - lambda_n  Thermal conductivities matrix in northface
%--------------------------------------------------------------------------


N = N1 + N2;        % Number of CVs in x direction
L = L1 + L2 + L3;   % Number of CVs in y direction

% Thermal conductivities on east face of external nodes
lambda_e = zeros(L+2,N+2);
lambda_e(:,1) = lambda(:,2);        % Column 1
lambda_e(:,end) = lambda(:,end);    % Not necessary, computed for completeness

% Thermal conductivities on north face of external nodes
lambda_n = zeros(L+2,N+2);
lambda_n(1,:) = lambda(2,:);        % Row 1
lambda_n(end,:) = lambda(end,:);    % Not necessary, computed for completeness

% Compute thermal conductivities of internal nodes using harmonic mean
for j = 2:L+1
    for i = 2:N+1
        lambda_e(j,i) = dx(i)/((x_cv(i)-x(i))/lambda(j,i) + (x(i+1)-x_cv(i))/lambda(j,i+1));
        lambda_n(j,i) = dy(j)/((y_cv(j)-y(j))/lambda(j,i) + (y(j+1)-y_cv(j))/lambda(j+1,i));
    end
end

% Thermal conductivities on west face
lambda_w = zeros(L+2,N+2);
lambda_w(:,1) = lambda(:,1);                % Not necessary, computed for completeness
lambda_w(:,2:end) = lambda_e(:,1:end-1);

% Thermal conductivities on south face
lambda_s = zeros(L+2,N+2);
lambda_s(1,:) = lambda(1,:);                % Not necessary, computed for completeness
lambda_s(2:end,:) = lambda_n(1:end-1,:);     

end