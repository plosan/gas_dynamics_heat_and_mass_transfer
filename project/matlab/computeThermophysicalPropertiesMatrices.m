function [rho, cp, lambda] = ...
    computeThermophysicalPropertiesMatrices(N1, N2, L1, L2, L3, rhoM, cpM, lambdaM)
%--------------------------------------------------------------------------
% Inputs: 
%   - N1        Number of CV in x direction for M1 and M3
%   - N2        Number of CV in x direction for M2 and M4
%   - L1        Number of CV in y direction for M1 and part of M2
%   - L2        Number of CV in y direction for part of M2 and part of M3
%   - rhoM      Densities of each material
%   - cpM       Specific heat of each material
%   - lambdaM   Thermal conductivity of each material
%--------------------------------------------------------------------------
% Outputs: 
%   - rho       (L1+L2+2)x(N1+N2+2) matrix containing the density of each node's material
%   - cp        (L1+L2+2)x(N1+N2+2) matrix containing the specific heat of each node's material
%   - lambda    (L1+L2+2)x(N1+N2+2) matrix containing the thermal conductivity of each node's material
%--------------------------------------------------------------------------

N = N1+N2;
L = L1+L2+L3;

rho = zeros(L+2,N+2);                   % Density matrix
rho(1:L1+1,1:N1+1) = rhoM(1);           % M1 density
rho(1:L1+L2+1,N1+2:end) = rhoM(2);      % M2 density
rho(L1+2:end,1:N1+1) = rhoM(3);         % M3 density
rho(L1+L2+2:end,N1+2:end) = rhoM(4);    % M4 density

cp = zeros(L+2,N+2);                % Specific heat matrix
cp(1:L1+1,1:N1+1) = cpM(1);         % M1 specific heat
cp(1:L1+L2+1,N1+2:end) = cpM(2);    % M2 specific heat
cp(L1+2:end,1:N1+1) = cpM(3);       % M3 specific heat
cp(L1+L2+2:end,N1+2:end) = cpM(4);  % M4 specific heat

lambda = zeros(L+2,N+2);                    % Thermal conductivity matrix
lambda(1:L1+1,1:N1+1) = lambdaM(1);         % M1 density
lambda(1:L1+L2+1,N1+2:end) = lambdaM(2);    % M2 density
lambda(L1+2:end,1:N1+1) = lambdaM(3);       % M3 density
lambda(L1+L2+2:end,N1+2:end) = lambdaM(4);  % M4 density

end