function outmatrix = vandermonde(x)

k = length(x);
vandermondematrix = zeros(k,k);
for i = 1:k
    for j = 1:k
        vandermondematrix(i,j) = x(i)^(k-j);
    end
end

outmatrix = vandermondematrix^(-1);

return
end