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
% 1.2.1 Number of control volumes
N1 = 25;            % Number of CV in x direction for M1 and M3
N2 = 6*N1/5;        % Number of CV in x direction for M2 and M4
L1 = 28;            % Number of CV in y direction for M1 and part of M2
L2 = 3*L1/4;        % Number of CV in y direction for part of M2 and part of M3
L3 = L1/4;          % Number of CV in y direction for part of M3 and M4
N = N1 + N2;        % Number of CV in x direction
L = L1 + L2 + L3;   % Number of CV in y direction
% 1.2.2 Time data
t0 = 0;                         % Initial time  [s] 
t_max = 5000;                    % Final time    [s]
t_step = 0.5;                     % Time step     [s]
n_max = 1+(t_max-t0)/t_step;    % Number of time steps
% 1.2.3 Convergencia criterion
delta = 1e-12;
itmax = 500;
% 1.2.4 Numerical integration
beta = 0.5; % Never 0
% filename = 'plots/verificacio/validacio_12.pdf';



%% 2. PREVIOUS COMPUTATIONS
% 2.1 Compute geometry (nodes position, distances, walls position, surfaces, volumes)
[delta_x1, delta_x2, delta_y1, delta_y2, delta_y3, ...
    x, y, x_cv, y_cv, dx, dy, Sx, Sy, V] = ...
    computeGeometry(p0, p1, p2, p3, W, N1, N2, L1, L2, L3, 0);

% 2.2 Compute thermophysical properties (density, specific heat, thermal conductivity)
[rho, cp, lambda] = computeThermophysicalPropertiesMatrices(N1, N2, L1, L2, L3, rhoM, cpM, lambdaM);

% 2.3 Compute thermal conductivity on faces using harmonic mean
[lambda_w, lambda_e, lambda_s, lambda_n] = computeThermalConductivities(N1, N2, L1, L2, L3, x, y, x_cv, y_cv, dx, dy, lambda);

% 2.4 Compute constant discretization coefficients
[A, b, d] = computeConstantDiscretizationCoefficients(N, L, ...
    t_step, beta, alpha_g, Tg, T_low, q_flow, dx, dy, Sx, Sy, V, rho, cp, lambda_w, lambda_e, lambda_s, lambda_n);

%% 3. INITIAL TEMPERATURE MAP
Tn = zeros((L+2)*(N+2),1)+T0;   % Vector of length (L+2)*(N+2) containing the temperature of each node at time n
Tn(1:N+2) = T_low;              % Initial condition

% Temperature matrix containing the temperature map for each time. Only
% used at the beginning of the development, to make sure the code worked
% properly
% T = zeros(L+2,N+2,n_max);   
% T(:,:,1) = T0;              % Initial temperature
% T(1,2:N+1,1) = T_low;       % Boundary condition on the lower wall

%% 4. COMPUTE NEXT TIME STEP AND VALIDATION
% Coordinates (j,i) of nodes to be validated. These coordinates strongly
% depend on the chosen spatial discretization. Make sure these nodes exist
% if you change the discretization.
nodes_val = [   % Coordinates (j,i) of nodes to be validated. 
    12  9;      % Internal node M1
    12  34;     % Internal node M2
    54  9;      % Internal node M3
    54  34;     % Internal node M4
    58  9;      % Upper wall node M3
    58  34;     % Upper wall node M4
    12  1;      % Left wall node M1
    54  1;      % Left wall node M3
];
nodes_val_internal = 1;     % First row of internal nodes in nodes_val
nodes_val_upper = 5;        % First row of upper wall nodes in nodes_val
nodes_val_left = 7;         % First row of left wall nodes in nodes_val
% Residues of validation. Each row corresponds to a time, each column 
% represents a node. Validation is performed at each iteration so that it 
% is not necessary to save a temperature map at each iteration. 
res = zeros(n_max-1, size(nodes_val,1));

for n = 1:n_max-1
    % COMPUTATION OF NON-CONSTANT DISCRETIZATION COEFFICIENTS
    % Internal nodes
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
    % Right wall nodes   
    b([2*(N+2):N+2:(L+1)*(N+2)]) = T_fun(n*t_step);    
    % SOLVE LINEAR SYSTEM MODIFIED-GAUSS-SEIDEL
    [T_MGS, it, res1] = modifiedGaussSeidel(A, b, Tn, N, L, itmax, delta);    
    % Singular nodes temperatures
    T_MGS(1) = T_MGS(2);
    T_MGS(N+2) = T_MGS(end-(N+2));
    T_MGS(1+(L+1)*(N+2)) = (T_MGS(2+(L+1)*(N+2)) + T_MGS(1+L*(N+2)))/2;
    T_MGS(end) = T_MGS((L+1)*(N+2));
    
    % VALIDATION    
    T1 = vectorToMatrix(Tn,L+2,N+2);        % Temperature at time n
    T2 = vectorToMatrix(T_MGS,L+2,N+2);     % Temperature at time n+1    
    % Validation of internal nodes
    for k = 1:nodes_val_upper-1        
        j = nodes_val(k,1);
        i = nodes_val(k,2);
        aW = lambda_w(j,i)*Sy(j)/dx(i-1);
        aE = lambda_e(j,i)*Sy(j)/dx(i);
        aS = lambda_s(j,i)*Sx(i)/dy(j-1);
        aN = lambda_n(j,i)*Sx(i)/dy(j);
        QP1 = aW*T1(j,i-1) + aE*T1(j,i+1) + aS*T1(j-1,i) + aN*T1(j+1,i) - (aW+aE+aS+aN)*T1(j,i);    % Net heat at time n
        QP2 = aW*T2(j,i-1) + aE*T2(j,i+1) + aS*T2(j-1,i) + aN*T2(j+1,i) - (aW+aE+aS+aN)*T2(j,i);    % Net heat at time n+1 n+1
        RHS = beta*QP2+(1-beta)*QP1;                            % Right hand side of the equation
        LHS = rho(j,i)*V(j,i)*cp(j,i)*(T2(j,i)-T1(j,i))/t_step; % Left hand side of the equation
        res(n,k) = abs(LHS-RHS);                                % Absolute value difference
    end
    % Validation of upper wall nodes
    for k = nodes_val_upper:nodes_val_left-1
        j = nodes_val(k,1);            
        i = nodes_val(k,2);
        res(n,k) = Sx(i)*abs(lambda_s(j,i)*(T2(j,i)-T2(j-1,i))/dy(j-1)-q_flow); % Absolute value difference        
    end
    % Validation of left wall nodes
    for k = nodes_val_left:size(nodes_val,1)
        j = nodes_val(k,1);            
        i = nodes_val(k,2);
        res(n,k) = Sy(j)*abs(alpha_g*(Tg-T2(j,i))+lambda_e(j,i)*(T2(j,i+1)-T2(j,i))/dx(i)); % Absolute value difference           
    end
    
    % UPDATE TEMPERATURE
    Tn = T_MGS;
end


%% 5. FINAL COMPUTATIONS AND PRINT RESULTS
h = figure(2);
hold on;
title(sprintf("\\textbf{Validaci\\'o amb esquema impl\\'icit, $\\Delta t = %.2f \\ \\mathrm{s}$}", t_step));
t_vec = [t0+t_step:t_step:t_max]';
plot(t_vec, res(:,1), 'r');
plot(t_vec, res(:,2), 'color', [255 149 117]/255);
plot(t_vec, res(:,3), 'g');
plot(t_vec, res(:,4), 'color', [0 119 4]/255);
plot(t_vec, res(:,5), 'b');
plot(t_vec, res(:,6), 'color', [115 147 255]/255);
plot(t_vec, res(:,7), 'm');
plot(t_vec, res(:,8), 'k');
xlabel("Temps $(\mathrm{s})$");
ylabel("Difer\`encia $| \mathrm{RHS}_k - \mathrm{LHS}_k | \ (\mathrm{W})$");
grid on;
box on;
set(gcf,'units','points','position',[100,100,500,375]);
legend_str = cell(8,1);
legend_str(1) = {sprintf("Node $(%d,%d)$", nodes_val(1,2), nodes_val(1,1))};
legend_str(2) = {sprintf("Node $(%d,%d)$", nodes_val(2,2), nodes_val(2,1))};
legend_str(3) = {sprintf("Node $(%d,%d)$", nodes_val(3,2), nodes_val(3,1))};
legend_str(4) = {sprintf("Node $(%d,%d)$", nodes_val(4,2), nodes_val(4,1))};
legend_str(5) = {sprintf("Node $(%d,%d)$", nodes_val(5,2), nodes_val(5,1))};
legend_str(6) = {sprintf("Node $(%d,%d)$", nodes_val(6,2), nodes_val(6,1))};
legend_str(7) = {sprintf("Node $(%d,%d)$", nodes_val(7,2), nodes_val(7,1))};
legend_str(8) = {sprintf("Node $(%d,%d)$", nodes_val(8,2), nodes_val(8,1))};
legend(legend_str);
hold off;

% Save plot as pdf
% set(h, 'Units', 'Centimeters');
% pos = get(h, 'Position');
% set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(h, filename, '-dpdf', '-r0');


