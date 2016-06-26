function outmatrix = powerseriescoeffs(x)
k = length(x);

vandermondematrix = zeros(k,k);
for i = 1:k
    for j = 1:k-1
        vandermondematrix(i,j) = x(i)^(k-j-1);
    end
    vandermondematrix(i,k) = x(i)^(-2);
end

outmatrix = vandermondematrix^(-1);

return
end