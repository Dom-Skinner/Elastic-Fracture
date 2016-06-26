function h_values_matrix = hprime_to_h_values(x)

n = length(x);
t = round(n/2);

h_values_matrix = zeros(n,2*n);

interpolate_matrix = linear_simpleinfty_interpolate(x);
interpolate_matrix = interpolate_matrix(2*n+1:4*n,n+1:2*n);

for j=1:n-1
    if j <= t-1
        for i=j+1:n
            h_values_matrix(i,j) = (2/3)*(x(i)^(3/2)-x(i-1)^(3/2));
            h_values_matrix(i,n+j) = 2*(x(i)^(1/2) - x(i-1)^(1/2));
        end
    else
        for i=j+1:n
            h_values_matrix(i,j) = (1/2)*(x(i)^(2)-x(i-1)^(2));
            h_values_matrix(i,n+j) = (x(i) - x(i-1));
        end
    end
end

h_values_matrix = h_values_matrix*interpolate_matrix;