clear;
close all;
clc;

%% 1. DATA INPUT
% 1.1 Physical data
% Geometrical data
R1 = 0.10;          % Cylinder internal radius  [m]
R2 = 0.20;          % Cylinder external radius  [m]
H = 1;              % Cylinder height           [m]

% Cylinder thermophysical data
qv = 1e4;           % Internal heat source      [W/m^3]
rho = 2700;         % Density                   [kg/m^3]
cp0 = 900;          % Specific heat             [J/kg K]
lambda0 = 50;       % Thermal conductivity      [W/m K]

% External fluids thermophysical data
alphaA0 = 1300;     % Heat transfer coefficient of fluid A  [W/m^2 K]
alphaB0 = 120;      % Heat transfer coefficient of fluid B  [W/m^2 K]

TA0 = 30;           % Average temperature of fluid A    [Celsius]
TB0 = 30;           % Average temperature of fluid B    [Celsius]

DeltaTA = 20;       % Amplitude of temperature oscillation of fluid A   [Celsius]
DeltaTB = 5;        % Amplitude of temperature oscillation of fluid A   [Celsius]

fA = 1/3600;        % Oscillation frequency for fluid A     [1/s]
fB = 1/(24*3600);   % Oscillation frequency for fluid B     [1/s]

TA = @(t) TA0 + DeltaTA*sin(2*pi*fA*t);    % Fluid A temperature   [Celsius]
TB = @(t) TB0 + DeltaTB*sin(2*pi*fB*t);    % Fluid B temperature   [Celsius]

lambda = @(T) 0*T+lambda0;
cp = @(T) 0*T+cp0;
alphaA = @(T) 0*T+alphaA0;
alphaB = @(T) 0*T+alphaB0;

% 1.2 Numerical data
N = 50;
delta = 1e-6;
beta = 0.5;
t_step = 1;
T0_star = 30;

%% 2. PREVIOUS COMPUTATIONS
% 2.1 Nodes position
DeltaR = (R2-R1)/N;
r_node = zeros(N+2,1);
r_node(1) = R1;
r_node(2:end-1) = [R1+DeltaR/2:DeltaR:R2-DeltaR/2];
r_node(end) = R2;

% 2.2 Distances
d = zeros(N+1,1); % d(i) is the distance between nodes i and i+1
for i = 1:N+1
    d(i) = r_node(i+1)-r_node(i);
end

% Surfaces
r_surf = [R1:DeltaR:R2]';
S = 2*pi*r_surf*H;

% Volumes
V = zeros(N+2,1);
for i = 2:N+1
    V(i) = pi*(r_surf(i)^2-r_surf(i-1)^2)*H;
end

%% 3. INITIAL MAP
t_max = 3600;    % Maximum time
M = t_max/t_step;   % Number of time instants
T = zeros(N+2,M+1);   % Temperature matrix
T(:,1) = T(:,1) + T0_star;  % Initial map
Q = zeros(N+2,M+1);         % Heat flows
Q(:,1) = 0;

%% 4. COMPUTE NEXT INSTANT OF TIME
n = 1;
while n < M+1
    % 4.1 ESTIMATION OF THE INITIAL TEMPERATURE MAP
    T_est = T(:,n);
    error = intmax;
    while error > delta
        % 4.2 COMPUTATION OF THE DISCRETIZATION COEFFICIENTS
        A = zeros(N+2,3);
        b = zeros(N+2,1);
        % Node 1
        A(1,1) = 0;
        A(1,3) = beta*lambda(T_est(2))*S(1)/d(1);
        A(1,2) = A(1,1) + A(1,3) + beta*alphaA(T_est(1))*S(1);
        QPn = alphaA(T(1,n))*(TA((n-1)*t_step)-T(1,n))*S(1)+...
            lambda(T(2,n))*(T(2,n)-T(1,n))*S(1)/d(1);
        b(1) = beta*alphaA(T_est(1))*TA(n*t_step)*S(1)+(1-beta)*QPn;
        % Internal nodes 2 <= i <= N+1
        for i = 2:N+1
            lambda_w = d(i-1)/((r_surf(i-1)-r_node(i-1))/lambda(T_est(i-1))...
                + (r_node(i)-r_surf(i-1))/lambda(T_est(i)));
            A(i,1) = beta*lambda_w*S(i-1)/d(i-1);
            lambda_e = d(i)/((r_surf(i)-r_node(i))/lambda(T_est(i)) + ...
                (r_node(i+1)-r_surf(i))/lambda(T_est(i+1)));
            A(i,3) = beta*lambda_e*S(i)/d(i);
            A(i,2) = A(i,1) + A(i,3) + rho*V(i)*cp(T_est(i))/t_step;

            lambda_w = d(i-1)/((r_surf(i-1)-r_node(i-1))/lambda(T(i-1,n)) +...
                (r_node(i)-r_surf(i-1))/lambda(T(i,n)));
            lambda_e = d(i)/((r_surf(i)-r_node(i))/lambda(T(i,n)) + ...
                (r_node(i+1)-r_surf(i))/lambda(T(i+1,n)));
            Qpn = -lambda_w*(T(i,n)-T(i-1,n))*S(i-1)/d(i-1) + ...
                lambda_e*(T(i+1,n)-T(i,n))*S(i)/d(i) + qv*V(i);
            b(i) = beta*qv*V(i) + rho*V(i)*cp(T(i,n))*T(i,n)/t_step + (1-beta)*QPn;
        end
        % Node N+2
        A(N+2,1) = beta*lambda(T_est(N+1))*S(N+1)/d(N+1);
        A(N+2,3) = 0;
        A(N+2,2) = A(N+2,1) + A(N+2,3) + beta*alphaB(T_est(N+2))*S(N+1);
        QPn = -lambda(T(N+1,n))*(T(N+2,n)-T(N+1,n))*S(N+1)/d(N+1) - ...
            alphaB(T(N+2,n))*(T(N+2,n)-TB((n-1)*t_step))*S(N+1);
        b(N+2) = beta*alphaB(T_est(N+2))*TB(n*t_step) + (1-beta)*QPn;
        % 4.3 SOLVE THE LINEAR SYSTEM USING TDMA
        T(:,n+1) = tdma(A, b);
        % 4.4 CHECK CONVERGENCE
        T_est = T(:,n+1);
        error = max(abs(T(:,n+1)-T_est));
    end
    % 
    % 5. NEXT INSTANT OF TIME
    n = n + 1;    
end

%% 5. PRINT RESULTS

% 5.1 Set interpreter latex
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% 5.2 3D Plot
[X,Y] = meshgrid(r_node, [0:t_step:t_max]);

figure(1);
hold on;
s = surf(X,Y,T');
s.EdgeColor = 'none';
xlabel("Posici\'o $(\mathrm{m})$");
ylabel("Temps $(\mathrm{s})$");
zlabel("Temperatura $(^\circ \mathrm{C})$");
set(gcf, 'units', 'pixels', 'position',[1000,400,1.2*560,1.2*420]);
ylim([0 3600]);
yticks([0:900:3600]);
box on;
grid on;
view([45 20]);
hold off;

% 5.3 T plot for different instants of time
t_plot = [0:900:3600]+1;
legend_str = cell(1,length(t_plot));
[asd, index] = max(max(T));
h3 = figure(3);
title("\textbf{Temperatura -- Posici\'o}");
hold on;
for i = 1:length(t_plot)
    legend_str(i) = {sprintf("$T(x)$ per $t = %d \\ \\mathrm{s}$", t_plot(i)-1)};
end
plot(r_node, T(:,t_plot(1)), 'r', 'LineWidth', 1);
plot(r_node, T(:,t_plot(2)), 'g', 'LineWidth', 1);
plot(r_node, T(:,t_plot(3)), 'b', 'LineWidth', 1);
plot(r_node, T(:,t_plot(4)), 'm', 'LineWidth', 1);
plot(r_node, T(:,t_plot(5)), 'c', 'LineWidth', 1);
xlabel("Posici\'o $(\mathrm{m})$");
ylabel("Temperatura $(^\circ \mathrm{C})$");
set(gcf, 'units', 'pixels', 'position',[1000,400,1.2*560,1.2*420]);
xlim([0.1 0.2]);
ylim([15 50]);
legend(legend_str);
box on;
grid on;
hold off;

% 5.4 Save plot as pdf
set(h3, 'Units', 'Centimeters');
pos = get(h3, 'Position');
set(h3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h3, 'figures/temperatura_posicio', '-dpdf', '-r0');























