function M = vectorToMatrix(v, n, m)

if length(v) ~= n*m
    error("Vector and matrix sizes do not match");
end

M = zeros(n,m);
for i = 1:n
    M(i,:) = v((i-1)*m+1:i*m);
end

end