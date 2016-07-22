%this function, for a given h', analytically integrates to give h
%use 1/sqrt(x) and linear representation for h'

%produces a matrix, which when multiplied by the n values of h',
%gives 3n coefficients for h.

function h_coefficient_matrix = hprime_to_h_s(x,s)

n = length(x);
t = round(n/2);

%gets the coefficients of g', h'
interpolate_matrix = linear_simpleinfty_interpolate_s(x,s);
interpolate_matrix = interpolate_matrix(2*n+1:4*n,n+1:2*n);

h_coefficient_matrix = zeros(3*n,2*n);


%gets the a and b coefficients - easy
for j=1:n
    if j <= t-1
        h_coefficient_matrix(3*(j-1)+1,j) = 1/(s+1);
        h_coefficient_matrix(3*(j-1)+2,n+j) = 1/s;
    else
        h_coefficient_matrix(3*(j-1)+1,j) = 1/2;
        h_coefficient_matrix(3*(j-1)+2,n+j) = 1;
    end
end

%gets the c coefficients.
for j=2:n
    for i = 1:j-1
        if i <= t-1
            h_coefficient_matrix(3*j,i) = (1/(s+1))*(x(i+1)^(s+1) - x(i)^(s+1));
            h_coefficient_matrix(3*j,n+i) = (1/s)*(x(i+1)^(s) - x(i)^(s));
        else
            h_coefficient_matrix(3*j,i) = (1/2)*(x(i+1)^2 - x(i)^2);
            h_coefficient_matrix(3*j,n+i) = x(i+1) - x(i);
        end
    end

    if j <= t-1
        h_coefficient_matrix(3*j,j) = -(1/(s+1))*x(j)^(s+1);
        h_coefficient_matrix(3*j,n+j) = -(1/s)*x(j)^(s);
    else
        h_coefficient_matrix(3*j,j) = -(1/2)*x(j)^2;
        h_coefficient_matrix(3*j,n+j) = -x(j);
    end
end

h_coefficient_matrix = h_coefficient_matrix*interpolate_matrix;