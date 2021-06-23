function [delta_x1, delta_x2, delta_y1, delta_y2, delta_y3, ...
    x, y, x_cv, y_cv, dx, dy, Sx, Sy, V] = ...
    computeGeometry(p0, p1, p2, p3, W, N1, N2, L1, L2, L3, plotGeometry)
%--------------------------------------------------------------------------
% Inputs: 
%   - p0            Point 0
%   - p1            Point 1
%   - p2            Point 2
%   - p3            Point 3
%   - W             Material depth
%   - N1            Number of CV in x direction for M1 and M3
%   - N2            Number of CV in x direction for M2 and M4
%   - L1            Number of CV in y direction for M1 and part of M2
%   - L2            Number of CV in y direction for part of M2 and part of M3
%   - L3            Number of CV in y direction for part of M3 and M4
%   - plotGeometry  If it is different from 0, plots the discretization
%--------------------------------------------------------------------------
% Outputs: 
%   - delta_x1      Width of CVs in M1 and M3
%   - delta_x2      Width of CVs in M2 and M4
%   - delta_y1      Height of CVs in M1 and part of M2
%   - delta_y2      Heigth of CVs in part of M2 and part of M3
%   - delta_y3      Height of CVs in part of M3 and M4
%   - x             Nodes position in x direction
%   - y             Nodes position in y direction
%   - x_cv          CV walls position in x direction
%   - y_cv          CV walls position in y direction
%   - dx            Distances between nodes in x direction. dx(i) is the distance from node i to i+1
%   - dy            Distances between nodes in y direction. dy(i) is the distance from node i to i+1
%   - Sx            Surfaces in x axis
%   - Sy            Surfaces in y axis
%   - V             Nodes volume matrix
%--------------------------------------------------------------------------

N = N1 + N2;
L = L1 + L2 + L3;

% 2.1 Geometry
delta_x1 = p1(1)/N1;            % Width of CVs in M1 and M3
delta_x2 = (p3(1)-p1(1))/N2;    % Width of CVs in M2 and M4
delta_y1 = p1(2)/L1;            % Height of CVs in M1 and part of M2
delta_y2 = (p2(2)-p1(2))/L2;    % Heigth of CVs in part of M2 and part of M3
delta_y3 = (p3(2)-p2(2))/L3;    % Height of CVs in part of M3 and M4
% 2.1.1 CV walls positions
x_cv = zeros(N+1,1);                            % CV walls position in x direction
x_cv(1:N1+1) = delta_x1*[0:1:N1];               % CV walls position in x direction, M1 and M3
x_cv(N1+2:end) = p1(1) + delta_x2*[1:1:N2];     % CV walls position in x direction, M2 and M4
y_cv = zeros(L+1,1);                            % CV walls position in y direction
y_cv(1:L1+1) = delta_y1*[0:1:L1];               % CV walls position in y direction, M1
y_cv(L1+2:L1+L2+1) = p1(2) + delta_y2*[1:1:L2]; % CV walls position in y direction, upper part of M2 and lower part of M3
y_cv(L1+L2+2:end) = p2(2) + delta_y3*[1:1:L3];  % CV walls position in y direction, upper part of M3 and M4
% 2.1.2 Nodes position
x = zeros(N+2,1);                                   % Nodes position in x direction
x(2:N1+1) = x_cv(2:N1+1)-delta_x1/2;                % Nodes position in x direction, M1 and M3
x(N1+2:end-1) = x_cv(N1+2:end) - delta_x2/2;        % Nodes position in x direction, M2 and M4
x(end) = x_cv(end);                                 % Last node position
y = zeros(L+2,1);                                   % Nodes position in y direction
y(2:L1+1) = y_cv(2:L1+1) - delta_y1/2;              % Nodes position in y direction, M1 and lower part of M2
y(L1+2:L1+L2+1) = y_cv(L1+2:L1+L2+1) - delta_y2/2;  % Nodes position in y direction, upper part of M2 and lower part of M3
y(L1+L2+2:end-1) = y_cv(L1+L2+2:end) - delta_y3/2;  % Nodes position in y direction, upper part of M3 and M4
y(end) = y_cv(end);                                 % Last node position

% 2.2 Distances
% 2.2.1 Distances in x direction
dx = zeros(N+1,1);                      % Distances between nodes in x direction. dx(i) is the distance from node i to i+1
dx(1) = delta_x1/2;                     % Distance from 1 to 2
dx(2:N1) = delta_x1;                    % Distances between nodes in x direction, M1 and M3
dx(N1+1) = (delta_x1+delta_x2)/2;       % Distance between the last node of M1 (M3) and the first of M2 (M4)
dx(N1+2:end-1) = delta_x2;              % Distance between nodes in x direction, M2 and M4
dx(end) = delta_x2/2;                   % Distance from node N+1 to N+2
% 2.2.2 Distances in y direction
dy = zeros(L+1,1);                      % Distances between nodes in y direction. dy(i) is the distance from node i to i+1 
dy(1) = delta_y1/2;                     % Distances from 1 to 2
dy(2:L1) = delta_y1;                    % Distances between nodes in y direction, M1 and lower part of M2
dy(L1+1) = (delta_y1+delta_y2)/2;       % Distance between the last node of M1 (lower part of M2) and the first node of the lower part of M3 (upper part of M2)
dy(L1+2:L1+L2) = delta_y2;              % Distances between nodes in y direction, lower part of M3 and upper part of M2
dy(L1+L2+1) = (delta_y2+delta_y3)/2;    % Distances between the last node of the lower part of M3 (upper part of M2) and the first node of the upper part of M3 (M4)
dy(L1+L2+2:end-1) = delta_y3;           % Distance between nodes in y direction, upper part of M3 and M4
dy(end) = delta_y3/2;                   % Distance from L+1 to L+2

% 2.3 Surfaces
% 2.3.1 Surfaces in x direction
Sx = zeros(N+2,1);                  % Surfaces along x axis (S and N)
Sx(2:N1+1) = delta_x1*W;            % Surfaces along x axis, M1 and M3
Sx(N1+2:end-1) = delta_x2*W;        % Surfaces along x axis, M2 and M4
% 2.3.1 Surfaces in y direction
Sy = zeros(L+2,1);                  % Surfaces along y axis (W and E)
Sy(2:L1+1) = delta_y1*W;            % Surfaces along y axis, M1 and lower part of M2
Sy(L1+2:L1+L2+1) = delta_y2*W;      % Surfaces along y axis, upper part of M2 and lower part of M3
Sy(L1+L2+2:end-1) = delta_y3*W;     % Surfaces along y axis, upper part of M3 and M4

% 2.4 Volumes
V = zeros(L+2,N+2);                                 % Nodes volume matrix
V(2:L1+1,2:N1+1) = delta_x1*delta_y1*W;             % Nodes in M1
V(L1+2:L1+L2+1,2:N1+1) = delta_x1*delta_y2*W;       % Nodes in lower part of M3
V(L1+L2+2:end-1,2:N1+1) = delta_x1*delta_y3*W;      % Nodes in upper part of M3
V(2:L1+1,N1+2:end-1) = delta_x2*delta_y1*W;         % Nodes in lower part of M2
V(L1+2:L1+L2+1,N1+2:end-1) = delta_x2*delta_y2*W;   % Nodes in upper part of M2
V(L1+L2+2:end-1,N1+2:end-1) = delta_x2*delta_y3*W;  % Nodes in M4

% 2.5 Geometry check plot
if plotGeometry ~= 0
    figure(1);
    hold on;
    rectangle('Position',[p0(1) p0(2) p3(1) p3(2)], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position',[p0(1) p0(2) p1(1) p1(2)], 'FaceColor', 'b', 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position',[p1(1) p0(2) p3(1)-p1(1) p2(2)], 'FaceColor', 'm', 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position',[p0(1) p1(2) p1(1) p3(2)-p1(2)], 'FaceColor', 'g', 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position',[p2(1) p2(2) p3(1)-p2(1) p3(2)-p2(2)], 'FaceColor', 'r', 'EdgeColor', 'k', 'LineWidth', 1);
    for i = 1:length(x_cv)
        plot([x_cv(i) x_cv(i)], [p0(2) p3(2)], '--k');
    end
    for i = 1:length(y_cv)
        plot([p0(1) p3(1)], [y_cv(i) y_cv(i)], '--k');
    end
    for i = 1:length(x)
        scatter(zeros(length(y),1)+x(i), y, 10, 'k', 'filled');
    end
    axis equal;
%     xlim([0 1.2]);
%     ylim([0 1]);
    hold off;
end

end


