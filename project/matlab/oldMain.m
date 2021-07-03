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
Tg = 33;        % Fluid temperature                 [ºC]
alpha_g = 9;    % Fluid heat transfer coefficient   [W/m^2 K]

% 1.1.4 Initial and boundary conditions
T0 = 8;                     % Initial temperature                           [ºC]
T_low = 23;                 % Lower wall node temperature                   [ºC]
Q_flow = 60;                % Upper wall heat flow per unit length of bar   [W/m]
q_flow = Q_flow/W;          % Upper wall heat flow per unit surface         [W/m^2]
T_fun = @(t) 8+0.005*t;     % Right wall temperature                        [ºC]

% 1.2 Numerical data
% 1.2.1 Number of control volumes
N1 = 1;             % Number of CV in x direction for M1 and M3
N2 = 1;             % Number of CV in x direction for M2 and M4
L1 = 1;             % Number of CV in y direction for M1 and part of M2
L2 = 1;             % Number of CV in y direction for part of M2 and part of M3
L3 = 1;             % Number of CV in y direction for part of M3 and M4

N1 = 20;
N2 = 6*N1/5;
L1 = N1;
L2 = 3*L1/5;
L3 = L1/5;

N = N1 + N2;        % Number of CV in x direction
L = L1 + L2 + L3;   % Number of CV in y direction
% 1.2.2 Time data
t0 = 0;                         % Initial time  [s] 
t_max = 100;                  % Final time    [s]
t_step = 1e-1;                     % Time step     [s]
n_max = 1+(t_max-t0)/t_step;    % Number of time steps
% 1.2.3 Convergencia criteria
delta = 1e-10;
itmax = 500;
% 1.2.4 Numerical integration
beta = 0.5;

%% 2. PREVIOUS COMPUTATIONS
% 2.1 Compute geometry (nodes position, distances, walls position, surfaces, volumes)
[delta_x1, delta_x2, delta_y1, delta_y2, delta_y3, ...
    x, y, x_cv, y_cv, dx, dy, Sx, Sy, V] = ...
    computeGeometry(p0, p1, p2, p3, W, N1, N2, L1, L2, L3, 1);

% 2.2 Compute thermophysical properties (density, specific heat, thermal conductivity)
[rho, cp, lambda] = computeThermophysicalPropertiesMatrices(N1, N2, L1, L2, L3, rhoM, cpM, lambdaM);

% 2.3 Compute thermal conductivity on faces using harmonic mean
[lambda_w, lambda_e, lambda_s, lambda_n] = computeThermalConductivities(N1, N2, L1, L2, L3, x, y, x_cv, y_cv, dx, dy, lambda);

% lambda_w = lambda;
% lambda_e = lambda;
% lambda_s = lambda;
% lambda_n = lambda;

%% 3. INITIAL TEMPERATURE MAP
T = zeros(L+2,N+2,n_max);
T(:,:,1) = T0;              % Initial temperature
T(1,2:N+1,1) = T_low;       % Boundary condition on the lower wall
T2 = T;

%% 3. COMPUTE DISCRETIZATION COEFFICIENTS
Tn = zeros((L+2)*(N+2),1)+T0;
Tn(1:N+2) = T_low;

A3 = zeros((L+2)*(N+2),5);
b3 = zeros((L+2)*(N+2),1);
d = zeros((L+2)*(N+2),1);

b3(2:N+1,1) = T_low;
b3(2+(L+1)*(N+2):(L+2)*(N+2)-1) = q_flow;
b3([N+3:N+2:1+L*(N+2)]) = alpha_g*Tg; 

% INTERNAL NODES COEFFICIENTS
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
        A3(id,1) = aS;
        A3(id,2) = aW;
        A3(id,3) = aP;
        A3(id,4) = aE;
        A3(id,5) = aN;
        d(id,1) = dP; 
    end
end

% LOWER AND UPPER WALL NODES
for i = 2:N+1
    % LOWER WALL NODES
    id = i;
    A3(id,3) = 1;
    % UPPER WALL NODES
    id = i+(L+1)*(N+2);
    aS = lambda_s(L+2,i)/dy(L+1);
    aP = aS;
    A3(id,1) = aS;
    A3(id,3) = aP;
end

% LEFT AND RIGHT WALL NODES
for j = 2:L+1
    % LEFT WALL NODES
    id = 1+(j-1)*(N+2);
    aE = lambda_e(j,1)/dx(1);
    aP = aE + alpha_g;
    A3(id,3) = aP;
    A3(id,4) = aE;
    % RIGHT WALL NODES
    id = j*(N+2);
    A3(id,3) = 1;
end


%% 4. COMPUTE NEXT TIME STEP

res = zeros(n_max,1);

for n = 1:n_max-1
%     Tn = matrixToVector(T(:,:,n));  % Temperature vector at time n
    % 4.1 COMPUTE DISCRETIZATION COEFFICIENTS
    A = zeros((N+2)*(L+2));
    A2 = zeros((N+2)*(L+2),5);
    b = zeros((N+2)*(L+2),1);
    % 4.1.1 INTERNAL NODES
    for j = 2:L+1
        for i = 2:N+1
            % Node P ID
            id = i+(j-1)*(N+2);
            % First version of the discretization coefficients
            aW = lambda_w(j,i)*Sy(j)/dx(i-1);
            aE = lambda_e(j,i)*Sy(j)/dx(i);
            aS = lambda_s(j,i)*Sx(i)/dy(j-1);
            aN = lambda_n(j,i)*Sx(i)/dy(j);
            % Net heat on node P at time n
            QPn = (aW*Tn(id-1)+aE*Tn(id+1)+aS*Tn(id-(N+2))+aN*Tn(id+(N+2))-(aW+aE+aS+aN)*Tn(id));
            % Discretization coefficients
            aW = beta*aW;
            aE = beta*aE;
            aS = beta*aS;
            aN = beta*aN;
            dP = rho(j,i)*V(j,i)*cp(j,i)/t_step;
            aP = dP+aW+aE+aS+aN;
            bP = dP*Tn(id)+(1-beta)*QPn;
            % Insert coefficients on matrix
            A(id,id) = aP;
            A(id,id-1) = -aW;
            A(id,id+1) = -aE;
            A(id,id-(N+2)) = -aS;
            A(id,id+(N+2)) = -aN;
            b(id,1) = bP;
            % Insert coefficients on matrix A2
            A2(id,1) = aS;
            A2(id,2) = aW;
            A2(id,3) = aP;
            A2(id,4) = aE;
            A2(id,5) = aN;
            % INDEPENDENT TERM            
            QPn = ((aW*Tn(id-1)+aE*Tn(id+1)+aS*Tn(id-(N+2))+aN*Tn(id+(N+2))-(aW+aE+aS+aN)*Tn(id)))/beta;
            b3(id) = d(id)*Tn(id)+(1-beta)*QPn;
        end
    end
    % 4.1.2 LOWER WALL AND UPPER WALL NODES
    for i = 2:N+1
        % LOWER WALL        
        % Insert coefficients on matrix
        A(i,i) = 1;
        b(i,1) = T_low;
        % Insert coefficients on matrix A2
        A2(i,3) = 1;
        
        % UPPER WALL
        % Node P ID
        id = i+(L+1)*(N+2);
        % Discretization coefficients
        aS = lambda_s(L+2,i)/dy(L+1);
        aP = aS;
        bP = q_flow;
        % Insert coefficients on matrix
        A(id,id) = aP;
        A(id,id-(N+2)) = -aS;
        b(id,1) = bP;
        % Insert coefficients on matrix A2
        A2(id,1) = aS;
        A2(id,3) = aP;
    end
    % 4.1.3 LEFT WALL AND RIGHT WALL NODES
    for j = 2:L+1
        % LEFT WALL
        % Node P ID
        id = 1+(j-1)*(N+2);
        % Discretization coefficients
        aE = lambda_e(j,1)/dx(1);
        aP = aE + alpha_g;
        bP = alpha_g*Tg;
        % Insert coefficients on matrix
        A(id,id) = aP;
        A(id,id+1) = -aE;
        b(id,1) = bP;
        % Insert coefficients on matrix A2
        A2(id,3) = aP;
        A2(id,4) = aE;
        
        % RIGHT WALL
        % Node P ID
        id = j*(N+2);
        % Insert coefficients on matrix
        A(id,id) = 1;
        b(id,1) = T_fun(n*t_step);
        % Insert coefficients on matrix A2
        A2(id,3) = 1;
        
        b3(id) = T_fun(n*t_step);
    end
    C = zeros((L+2)*(N+2),(L+2)*(N+2)+1);
    C(:,1:end-1) = A;
    C(:,end) = b;
    
    
    
    % 4.4 SOLVE LINEAR SYSTEM Mod-GS
    [T_MGS, it, res1] = modifiedGaussSeidel(A3, b3, Tn, N, L, itmax, delta);
    
    
    % 4.4.1 SINGULAR NODES TEMPERATURES
    T_MGS(1) = T_MGS(2);
    T_MGS(N+2) = T_MGS(end-(N+2));
    T_MGS(1+(L+1)*(N+2)) = (T_MGS(2+(L+1)*(N+2)) + T_MGS(1+L*(N+2)))/2;
    T_MGS(end) = T_MGS((L+1)*(N+2));
%     fprintf("%10s = %.5e\n", "T_MGS(end)", T_MGS(end));
    
%     res4 = norm(b-b3);
%     fprintf("%10s = %.5e\n", "res4", res4);
    
    % 4.2 ERASE SINGULAR NODES
    singularNodes = [1,N+2,1+(L+1)*(N+2),(L+2)*(N+2)];
    A(singularNodes,:) = [];
    A(:,singularNodes) = [];
    b(singularNodes,:) = [];
%     Tn(singularNodes,:) = [];
    A2(singularNodes,:) = [];
    % 4.3 SOLVE LINEAR SYSTEM LU
    T_new = linsolve(A,b);
    
    
    % COMPUTE RESIDUE
    
    
    % 4.4 COMPUTE SINGULAR NODES TEMPERATURE
    
    T_new2 = [0; T_new(1:N,1); 0; T_new(N+1:end-N,1); 0; T_new(end-N+1:end,1); 0];
    T_new2(1) = T_new(1);
    T_new2(N+2) = T_new(2*N+2);
    T_new2(1+(L+1)*(N+2)) = (T_new(end-N+1)+T_new(end-(N+N+2)+1))/2;
    T_new2(end) = T_new(end-N);    
    T_new = T_new2;
    
%     [res2,index] = max(abs(T_MGS-T_new));
    
    
    res(n+1) = max(abs(T_MGS-T_new));    
    
    % 4.5 SAVE TEMPERATURE
    T(:,:,n+1) = vectorToMatrix(T_new,L+2,N+2); 
    T2(:,:,n+1) = vectorToMatrix(T_MGS,L+2,N+2); 
    % 4.6 UPDATE TEMPERATURE
    Tn = T_new;
%     break;
    fprintf("%5d%10s%5.2f%10s%10.5e\n", n, "", n*t_step, "", res(n+1));
end

K = T(:,:,end);

%% 7. VERIFICATION
% res = zeros(L+2,N+2,n_max);
% 
% for n = 2:n_max-1
%     % 7.1 INTERNAL NODES
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
%             res(j,i,n) = abs(LHS-RHS);
%         end
%     end
%     % 7.2 UPPER WALL NODES
%     for i = 2:N+1
%         j = L+2;
%         res(j,i) = abs(q_flow-lambda_s(j,i)*(T(j,i,n)-T(j-1,i,n))/dy(j-1));
%     end
%     % 7.3 LEFT WALL NODES
%     for j = 2:L+1
%         i = 1;
%         res(j,i) = alpha_g*(Tg-T(j,i,n))+lambda_e(j,i)*(T(j,i+1,n)-T(j,i,n))/dx(i);
%     end
% end



%% 8. PLOT
% x_plot = zeros((N+2)*(L+2),1);
% y_plot = x_plot;
% for j = 1:L+2
%     x_plot((j-1)*(N+2)+1:j*(N+2),1) = x;
%     y_plot((j-1)*(N+2)+1:j*(N+2),1) = y(j);
% end
% T(1,1,end) = (T(2,1,end)+T(1,2,end)+T(2,2,end))/3;
% T(1,N+2,end) = (T(1,N+1,end)+T(2,N+2,end)+T(2,N+1,end))/3;
% T(L+2,1,end) = (T(L+1,1,end)+T(L+2,2,end)+T(L+1,2,end))/3;
% T(L+2,N+2,end) = (T(L+2,N+1,end)+T(L+1,N+2,end)+T(L+1,N+1,end))/3;





[X,Y] = meshgrid(x,y);
figure(2);
hold on;
s = surf(X,Y,T(:,:,end));
set(s, 'FaceColor', 'interp');
set(s, 'EdgeColor', 'none');
hold off;


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

