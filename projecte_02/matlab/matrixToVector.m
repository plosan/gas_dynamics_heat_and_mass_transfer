function v = matrixToVector(A)

L = size(A,1);      % Number of rows
N = size(A,2);      % Number of columns
v = zeros(L*N,1);   % Output vector
for j = 1:L
    v((j-1)*N+1:j*N) = A(j,:)';
end

end