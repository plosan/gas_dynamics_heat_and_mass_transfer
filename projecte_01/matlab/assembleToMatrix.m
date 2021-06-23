function A = assembleToMatrix(T)
%--------------------------------------------------------------------------
% Description: this function creates a tri-diagonal n x n matrix from a 
% [n x 3] matrix 
%--------------------------------------------------------------------------
% Inputs: 
%   - T         Matrix containing the diagonals [n x 3]
%--------------------------------------------------------------------------
% Outputs: 
%   - A         Tri-digonal matrix. The main diagonal is T(:,2), the lower
%               diagonal is T(:,1) and the upper diagonal is T(:,2)
%--------------------------------------------------------------------------

n = size(T,1);      % Matrix dimension
A = zeros(n,n);     % Output matrix

T(:,1) = -T(:,1);
T(:,3) = -T(:,3);

% Create tri-diagonal matrix
% Fill rows 2 to n-1
for i = 2:n-1
    A(i,i-1:i+1) = T(i,:);
end
% Fill rows 1 and n
for j = 1:2
    A(1,j) = T(1,j+1);
    A(n,n-2+j) = T(n,j);
end

end