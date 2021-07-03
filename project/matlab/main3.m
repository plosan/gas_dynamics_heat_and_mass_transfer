clear;
close all;
clc;
format long;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1. DATA INPUT
% 1.1 Physical data
% 1.1.1 Geometry
p0 = [0.00 0.00];   % Point 0       [m]
p1 = [0.50 0.40];   % Point 1       [m]
p2 = [0.50 0.70];   % Point 2       [m]
p3 = [1.10 0.80];   % Point 3       [m]
W = 1;              % Bar length    [m]
% 1.1.2 Materials thermophysical properties (M1 to M4)
rhoM = [1500 1600 1900 2500];   % Density               [kg/m^3] 
cpM = [750 770 810 930];        % Specific heat         [J/kg K]
lambdaM = [170 140 200 140];    % Thermal conductivity  [W/m K]

% rhoM = [1500 1500 1500 1500];   % Density               [kg/m^3] 
% cpM = [750 750 750 750];        % Specific heat         [J/kg K]
% lambdaM = [170 170 170 170];    % Thermal conductivity  [W/m K]
% 1.1.3 Fluid thermophysical properties
alpha_g = 9;    % Fluid heat transfer coefficient   [W/m^2 K]
Tg = 33;        % Fluid temperature                 [ºC]

% 1.1.4 Initial and boundary conditions
T0 = 8;                     % Initial temperature                           [ºC]
T_low = 23;                 % Lower wall node temperature                   [ºC]
Q_flow = 60;                % Upper wall heat flow per unit length of bar   [W/m]
q_flow = Q_flow/W;          % Upper wall heat flow per unit surface         [W/m^2]
T_fun = @(t) 8+0.005*t;     % Right wall temperature                        [ºC]

% 1.2 Numerical data
% 1.2.1 Control surfaces position
x_cv = [0:0.05:0.50 0.505:0.01:1.09 1.1]';
y_cv = [0:0.05:0.40 0.405:0.01:0.70 0.70 0.705:0.01:0.80 0.80]';

% 1.2.2 Time data
t0 = 0;                         % Initial time  [s] 
t_max = 10000;                    % Final time    [s]
t_step = 0.5;                     % Time step     [s]
n_max = 1+(t_max-t0)/t_step;    % Number of time steps
% 1.2.3 Convergencia criterion
delta = 1e-11;
itmax = 500;
% 1.2.4 Numerical integration
beta = 0.5; % Never 0

%% 2. PREVIOUS COMPUTATIONS
% 2.1 Compute geometry (nodes position, distances, walls position, surfaces, volumes)
[N, N1, N2, L, L1, L2, L3, x, y, dx, dy, Sx, Sy, V] = computeGeometry2(p0, p1, p2, p3, x_cv, y_cv, W, 1);

% 2.2 Compute thermophysical properties (density, specific heat, thermal conductivity)
[rho, cp, lambda] = computeThermophysicalPropertiesMatrices(N1, N2, L1, L2, L3, rhoM, cpM, lambdaM);
% 
% 2.3 Compute thermal conductivity on faces using harmonic mean
[lambda_w, lambda_e, lambda_s, lambda_n] = computeThermalConductivities(N1, N2, L1, L2, L3, x, y, x_cv, y_cv, dx, dy, lambda);

% 2.4 Compute constant discretization coefficients
[A, b, d] = computeConstantDiscretizationCoefficients(N, L, ...
    t_step, beta, alpha_g, Tg, T_low, q_flow, dx, dy, Sx, Sy, V, rho, cp, lambda_w, lambda_e, lambda_s, lambda_n);

% 2.5 Search nodes of interest indexs
i_node1 = find(x == 0.65);  % i index of node with x = 0.65
j_node1 = find(y == 0.56);  % j index of node with y = 0.56
i_node2 = find(x == 0.74);  % i index of node with x = 0.74
j_node2 = find(y == 0.72);  % j index of node with y = 0.72

%% 3. INITIAL TEMPERATURE MAP
Tn = zeros((L+2)*(N+2),1)+T0;   % Temperature vector
Tn(1:N+2) = T_low;              % Initial condition

T_nodes = zeros(n_max,3);       % Matrix containing the temperatures of nodes of interest
T_nodes(1,:) = [0 T0 T0];       % Initial temperature 

%% 4. COMPUTE NEXT TIME STEP
for n = 1:n_max-1
    % COMPUTATION OF NON-CONSTANT DISCRETIZATION COEFFICIENTS
    % INTERNAL NODES
    for j = 2:L+1
        for i = 2:N+1
            id = i+(j-1)*(N+2);     % Node id
            aS = A(id,1);           % aS discretization coefficient
            aW = A(id,2);           % aW discretization coefficient
            aP = A(id,3);           % aP discretization coefficient
            aE = A(id,4);           % aE discretization coefficient
            aN = A(id,5);           % aN discretization coefficient
            % Net heat flow at time n
            QPn = ((aW*Tn(id-1)+aE*Tn(id+1)+aS*Tn(id-(N+2))+aN*Tn(id+(N+2))-(aW+aE+aS+aN)*Tn(id)))/beta;
            % Independent term
            b(id) = d(id)*Tn(id)+(1-beta)*QPn;
        end
    end    
    % RIGHT WALL NODES    
    b([2*(N+2):N+2:(L+1)*(N+2)]) = T_fun(n*t_step);    
    % SOLVE LINEAR SYSTEM Mod-GS
    [T_MGS, it, res1] = modifiedGaussSeidel(A, b, Tn, N, L, itmax, delta);    
    % SINGULAR NODES TEMPERATURES
    T_MGS(1) = T_MGS(2);
    T_MGS(N+2) = T_MGS(end-(N+2));
    T_MGS(1+(L+1)*(N+2)) = (T_MGS(2+(L+1)*(N+2)) + T_MGS(1+L*(N+2)))/2;
    T_MGS(end) = T_MGS((L+1)*(N+2));
    % SAVE NODES OF INTEREST TEMPERATURE
    id_node1 = i_node1 + (j_node1-1)*(N+2);
    id_node2 = i_node2 + (j_node2-1)*(N+2);
    T_nodes(n+1,1) = n*t_step;
    T_nodes(n+1,2) = T_MGS(id_node1);
    T_nodes(n+1,3) = T_MGS(id_node2); 
    % 4.6 UPDATE TEMPERATURE
    Tn = T_MGS;
    fprintf("%5d%10s%.5f\n", n, "", n*t_step);
end

%% 5. FINAL COMPUTATIONS AND PRINT RESULTS
fileID = fopen('LOPEZ.txt','w');
fprintf(fileID,"%10d%5s%20.15f%5s%20.15f\n",0, "", T0, "", T0);
for i = 3:2:size(T_nodes,1)
    fprintf(fileID, "%10d%5s%20.15f%5s%20.15f\n", T_nodes(i,1), "", T_nodes(i,2), "", T_nodes(i,3));
end
fclose(fileID);

