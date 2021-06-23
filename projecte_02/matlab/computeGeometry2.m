function [N, N1, N2, L, L1, L2, L3, x, y, dx, dy, Sx, Sy, V] = computeGeometry2(p0, p1, p2, p3, x_cv, y_cv, W, plotGeometry)
%--------------------------------------------------------------------------
% Inputs: 
%   - p0            Point 0
%   - p1            Point 1
%   - p2            Point 2
%   - p3            Point 3
%   - x_cv          CV walls position in x direction
%   - y_cv          CV walls position in y direction
%   - W             Bar length
%   - plotGeometry  If it is different from 0, plots the discretization
%--------------------------------------------------------------------------
% Outputs: 
%   - N             Number of CV in x direction
%   - N1            Number of CV in x direction for M1 and M3
%   - N2            Number of CV in x direction for M2 and M4
%   - L             Number of CV in y direction
%   - L1            Number of CV in y direction for M1 and part of M2
%   - L2            Number of CV in y direction for part of M2 and part of M3
%   - L3            Number of CV in y direction for part of M3 and M4
%   - x             Nodes position in x direction
%   - y             Nodes position in y direction
%   - dx            Distances between nodes in x direction. dx(i) is the distance from node i to i+1
%   - dy            Distances between nodes in y direction. dy(i) is the distance from node i to i+1
%   - Sx            Surfaces in x axis
%   - Sy            Surfaces in y axis
%   - V             Nodes volume matrix
%   - plotGeometry  If it is different from 0, plots the discretization
%--------------------------------------------------------------------------

% Compute nodes position
x = zeros(length(x_cv)+1,1);
y = zeros(length(y_cv)+1,1);

x(1) = x_cv(1);
x(end) = x_cv(end);
for i = 2:length(x_cv)
    x(i) = (x_cv(i-1)+x_cv(i))/2;
end

y(1) = y_cv(1);
y(end) = y_cv(end);
for j = 2:length(y_cv)
    y(j) = (y_cv(j-1)+y_cv(j))/2;
end

% Compute distances
dx = zeros(length(x)-1,1);
for i = 1:length(dx)
    dx(i) = x(i+1)-x(i);
end
dy = zeros(length(y)-1,1);
for j = 1:length(dy)
    dy(j) = y(j+1)-y(j);
end

% Compute number of nodes in x direction
N = length(x)-2;
i = 1;
while(x_cv(i) <= 0.50)
    i = i + 1;
end
N1 = i-2;
N2 = N - N1;

% Compute number of nodes in y direction
L = length(y)-2;
j = 1;
while(y_cv(j) <= 0.40)
    j = j + 1;
end
L1 = j - 2;
while(y_cv(j) <= 0.70)
    j = j + 1;
end
L2 = j - L1 - 2;
L3 = L - L1 - L2;

% Compute surfaces
Sx = dx*W;
Sy = dy*W;

% Compute volumnes
V = zeros(L+2,N+2);
for j = 2:L+1
    for i = 2:N+1
        V(j,i) = W*(x_cv(i)-x_cv(i-1))*(y_cv(j)-y_cv(j-1));
    end
end

% Geometry plot for testing purposes
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
    hold off;
end

end


