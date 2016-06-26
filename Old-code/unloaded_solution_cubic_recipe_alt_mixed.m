%this gives a solution using the mixed cubic spline.

n = 250;
t = round(n/2);
endpoint = 30;

Pload = 0;
Mload = 1;

%data points
x = tan((1:n)*(atan(endpoint))/n);
%points to measure pressure at
z = tan((0.5:1:n+0.5)*(atan(endpoint))/n);

%first we need to find a matrix that interpolates any y1, y2, ...,
%yn, y1'', yn'' and gives us our desired coefficients

%we use DIFFERENT representations for g and h (linear, quadratic + error
%terms respectively) and then put boundary conditions on the form of the
%LAST panel

%creates our coefficient matrix
g_interpol_matrix = zeros(4*(n-1),n+2);
for i=1:n+2
    if i == n+1
        y = zeros(1,n);
        y1ddash = 1;
        ynddash = 0;
        g_interpol_matrix(1:4*(n-1),n+1) = cubic_sqrt_spline_recipe_mixed_x1(x,y,y1ddash,ynddash);
    elseif i == n+2
        y = zeros(1,n);
        y1ddash = 0;
        ynddash = 1;
        g_interpol_matrix(1:4*(n-1),n+2) = cubic_sqrt_spline_recipe_mixed_x1(x,y,y1ddash,ynddash);
    else
        y = zeros(1,n);
        y(i) = 1;
        y1ddash = 0;
        ynddash = 0;
        g_interpol_matrix(1:4*(n-1),i) = cubic_sqrt_spline_recipe_mixed_x1(x,y,y1ddash,ynddash);
    end
end


h_interpol_matrix = zeros(4*(n-1),n+2);
for i=1:n+2
    if i == n+1
        y = zeros(1,n);
        y1ddash = 1;
        ynddash = 0;
        h_interpol_matrix(1:4*(n-1),n+1) = cubic_sqrt_spline_recipe_mixed_x2(x,y,y1ddash,ynddash);
    elseif i == n+2
        y = zeros(1,n);
        y1ddash = 0;
        ynddash = 1;
        h_interpol_matrix(1:4*(n-1),n+2) = cubic_sqrt_spline_recipe_mixed_x2(x,y,y1ddash,ynddash);
    else
        y = zeros(1,n);
        y(i) = 1;
        y1ddash = 0;
        ynddash = 0;
        h_interpol_matrix(1:4*(n-1),i) = cubic_sqrt_spline_recipe_mixed_x2(x,y,y1ddash,ynddash);
    end
end

coefficient_matrix = zeros(8*(n-1),2*n+4);
%the first n+2 columns are dedicated to g
%second n+2 columns dedicated to h
coefficient_matrix(1:4*(n-1),1:n+2) = g_interpol_matrix;
coefficient_matrix(4*(n-1)+1:8*(n-1),n+3:2*n+4) = h_interpol_matrix;

%h_infty_coeffs(1,n) = 1;
%h_infty_coeffs(2,1:n+2) = [1, 0, -x(n)^(-2), -2*x(n)^(-3)]*interpol_matrix(4*(n-2)+1:4*(n-2)+4,:);
%h_infty_coeffs(3,1:n+2) = [0, 0, 2*x(n)^(-3), 6*x(n)^(-4)]*interpol_matrix(4*(n-2)+1:4*(n-2)+4,:);
%h_infty_coeffs(4,n-1) = 1;
%h_infty_coeffs(5,1:n+2) = [1, 0, -x(n-1)^(-2), -2*x(n-1)^(-3)]*interpol_matrix(4*(n-3)+1:4*(n-3)+4,:);
%h_infty_coeffs = [x(n)^2, x(n), 1, x(n)^(-1), x(n)^(-2); 2*x(n), 1, 0, -x(n)^(-2), -2*x(n)^(-3); ...
%    2, 0, 0, 2*x(n)^(-3), 6*x(n)^(-4); x(n-1)^2, x(n-1), 1, x(n-1)^(-1), x(n-1)^(-2); ...
%    2*x(n-1), 1, 0, -x(n-1)^(-2), -2*x(n-1)^(-3)]\h_infty_coeffs;

%creates the matrix!!
%creates the matrix to measure all the pressures
kernel_matrix = zeros(2*n+2,8*(n-1));
for j = 1:n+1
    [pressurerow, shearrow] = pressure_cubic_recipe_alt_mixed(z(j),x);
    kernel_matrix(j, :) = pressurerow;
    kernel_matrix(n+1+j,:) = shearrow;
end

%sets conditions at infinity: in order
%x coeff of g, 1 coeff of g
%x^2 coeff of h, x coeff of h
g_infty_conditions = g_interpol_matrix(4*(n-1)-2,:);
h_infty_conditions = h_interpol_matrix(4*(n-1)-2,:);

equation_matrix = zeros(2*n+4,2*n+4);
equation_matrix(1:2*n+2,1:2*n+4) = kernel_matrix*coefficient_matrix;
equation_matrix(2*n+3,1:n+2) = g_infty_conditions;
equation_matrix(2*n+4,n+3:2*n+4) = h_infty_conditions;

inhomog = zeros(2*n+4,1);
%inhomog(2*n+1) = 0;              %x coeff of g
inhomog(2*n+3) = (-Pload+6*Mload);       %1 coeff of g
%inhomog(2*n+3) = 0;               %x^2 coeff of h
inhomog(2*n+4) = (12*Mload);        %x coeff of h

%normalizes it / dimensionalizes
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

%finds the solution!!!
soln_vector = equation_matrix\inhomog;
conditioning = rcond(equation_matrix);

%finds KI, KII
KI = h_interpol_matrix(4,1:n+2)*soln_vector(n+3:2*n+4);
KII = g_interpol_matrix(4,1:n+2)*soln_vector(1:n+2);

gprime = soln_vector(1:n);
for i = 1:t
    gprime(i) = (1/sqrt(x(i)))*gprime(i);
end
for i = t+1:n
    gprime(i) = (1/x(i)^2)*gprime(i);
end

hprime = soln_vector(n+3:2*n+2);
for i = 1:t
    hprime(i) = (1/sqrt(x(i)))*hprime(i);
end
for i = t+1:n
    hprime(i) = (1/x(i))*hprime(i);
end

%%measure = tan((0.05:0.1:n-0.05)*atan(endpoint)/n);
%m = length(measure);
%pressure = zeros(1,m);
%shear = zeros(1,m);

%for i = 1:m
%    [pressurerow,shearrow] = pressure_cubic_recipe_alt_mixed(measure(i),x);
%    pressure(i) = pressurerow*coefficient_matrix*soln_vector;
%    shear(i) = shearrow*coefficient_matrix*soln_vector;
%end
        
        
        
        
        