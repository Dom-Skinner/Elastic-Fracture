%PROBLEM: given y'(0), y(0), y(1), y(2) ..., y(k), find a quadratic spline.
%points specified at x(0), x(1), ..., x(k)
%want a matrix to do this automatically.

%want to output 3k numbers.

function quadinterpol_matrix = quadratic_interpolation(x)

k = length(x)-1;
quadinterpol_matrix = zeros(3*k, k+2);

function coeffs = given_quadratic_interpolation(y, x, k)

coeffs = zeros(1,3*k);

coeffs([1,k+1,2*k+1]) = [x(1), 1, x(1)^2; 1, 0, 2*x(1); x(2), 1, x(2)^2]^(-1) ...
    * [y(2); y(1); y(3)];

for i = 2:k
    coeffs([i,k+i,2*k+i]) = [x(i),1,x(i)^2; 1, 0, 2*x(i); x(i+1),1,x(i+1)^2]^(-1) ...
        * [y(i+1); 2*coeffs(2*k+i-1)*x(i)+coeffs(i-1); y(i+2)];
end

end

for j = 1:k+2
    v = zeros(1,k+2);
    v(j) = 1;
    matrix_coeffs = given_quadratic_interpolation(v, x, k);
    
    quadinterpol_matrix(:,j) = matrix_coeffs;
end

return
end
