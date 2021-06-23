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
L1 = 4*N1/5;        % Number of CV in y direction for M1 and part of M2
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
delta = 1e-11;
itmax = 500;
% 1.2.4 Numerical integration
beta = 0.5; % Never 0

%% 2. PREVIOUS COMPUTATIONS
% 2.1 Compute geometry (nodes position, distances, walls position, surfaces, volumes)
[delta_x1, delta_x2, delta_y1, delta_y2, delta_y3, ...
    x, y, x_cv, y_cv, dx, dy, Sx, Sy, V] = ...
    computeGeometry(p0, p1, p2, p3, W, N1, N2, L1, L2, L3, 1);

% 2.2 Compute thermophysical properties (density, specific heat, thermal conductivity)
[rho, cp, lambda] = computeThermophysicalPropertiesMatrices(N1, N2, L1, L2, L3, rhoM, cpM, lambdaM);

% 2.3 Compute thermal conductivity on faces using harmonic mean
[lambda_w, lambda_e, lambda_s, lambda_n] = computeThermalConductivities(N1, N2, L1, L2, L3, x, y, x_cv, y_cv, dx, dy, lambda);

% 2.4 Compute constant discretization coefficients
[A, b, d] = computeConstantDiscretizationCoefficients(N, L, ...
    t_step, beta, alpha_g, Tg, T_low, q_flow, dx, dy, Sx, Sy, V, rho, cp, lambda_w, lambda_e, lambda_s, lambda_n);

%% 3. INITIAL TEMPERATURE MAP
T = zeros(L+2,N+2,n_max);
T(:,:,1) = T0;              % Initial temperature
T(1,2:N+1,1) = T_low;       % Boundary condition on the lower wall

%% 3. INITIAL TEMPERATURE MAP
Tn = zeros((L+2)*(N+2),1)+T0;
Tn(1:N+2) = T_low;

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
    % 4.5 SAVE TEMPERATURE
    T(:,:,n+1) = vectorToMatrix(T_MGS,L+2,N+2);    
    % 4.6 UPDATE TEMPERATURE
    Tn = T_MGS;
    fprintf("%5d%10s%.5f\n", n, "", n*t_step);
end


%% 7. VERIFICATION
% res = zeros(n_max,1);
% for n = 2:n_max-1
%     res_n = 0;
%     % INTERNAL NODES
%     for j = 2:L+1
%         for i = 2:N+1
%             % Coefficients
%             aW = lambda_w(j,i)*Sy(j)/dx(i-1);
%             aE = lambda_e(j,i)*Sy(j)/dx(i);
%             aS = lambda_s(j,i)*Sx(i)/dy(j-1);
%             aN = lambda_n(j,i)*Sx(i)/dy(j);
%             QP1 = aW*T(j,i-1,n)+aE*T(j,i+1,n)+aS*T(j-1,i,n)+aN*T(j+1,i,n)-(aW+aE+aS+aN)*T(j,i,n);
%             QP2 = aW*T(j,i-1,n+1)+aE*T(j,i+1,n+1)+aS*T(j-1,i,n+1)+aN*T(j+1,i,n+1)-(aW+aE+aS+aN)*T(j,i,n+1);
%             LHS = rho(j,i)*V(j,i)*cp(j,i)*(T(j,i,n+1)-T(j,i,n))/t_step;
%             RHS = beta*QP2+(1-beta)*QP1;
%             res_i = abs(LHS-RHS);
%             res_n = max(res_i,res_n);
%         end
%     end
%     % UPPER WALL NODES
%     for i = 2:N+1
%         j = L+2;
%         res_u = abs(q_flow-lambda_s(j,i)*(T(j,i,n)-T(j-1,i,n))/dy(j-1));
%         res_n = max(res_u,res_n);        
%     end
%     % LEFT WALL NODES
%     for j = 2:L+1
%         i = 1;
%         res_l = abs(alpha_g*(Tg-T(j,i,n))+lambda_e(j,i)*(T(j,i+1,n)-T(j,i,n))/dx(i));
%         res_n = max(res_l,res_n);
%     end
%     % SAVE MAX DIFFERENCE IN VECTOR
%     res(n) = res_n;
% end


%% 8. PLOT


T_end = matrixToVector(T(:,:,end));
x_plot = zeros((N+2)*(L+2),1);
y_plot = x_plot;
for j = 1:L+2
    x_plot((j-1)*(N+2)+1:j*(N+2),1) = x;
    y_plot((j-1)*(N+2)+1:j*(N+2),1) = y(j);
end

nx = length(unique(x_plot)); 
ny = length(unique(y_plot)); 
X = reshape(x_plot,nx,ny); 
Y = reshape(y_plot,nx,ny); 
Z = reshape(T_end,nx,ny); 

figure(3);
hold on;
s = pcolor(X,Y,Z);
% set(s, 'FaceColor', 'interp');
set(s, 'EdgeColor', 'none');
hold off;



%% 9. FINAL COMPUTATIONS AND PRINT RESULTS

