%this gives a solution using the mixed cubic spline.

n = 100;
t = round(n/2);
endpoint = 30;

Pload = 0;
Mload = 1;

%data points
x = tan((0:n-1)*(atan(endpoint))/(n-1));
%points to measure pressure at
z = tan((0.5:1:n-0.5)*(atan(endpoint))/(n-1));

%first we need to find a matrix that interpolates any y1, y2, ...,
%yn, y1'', yn'' and gives us our desired coefficients

%creates our coefficient matrix
interpol_matrix = zeros(4*(n-1),n+2);
for i=1:n+2
    if i == n+1
        y = zeros(1,n);
        y1ddash = 1;
        ynddash = 0;
        interpol_matrix(1:4*(n-1),n+1) = cubic_sqrt_spline_recipe(x,y,y1ddash,ynddash);
    elseif i == n+2
        y = zeros(1,n);
        y1ddash = 0;
        ynddash = 1;
        interpol_matrix(1:4*(n-1),n+2) = cubic_sqrt_spline_recipe(x,y,y1ddash,ynddash);
    else
        y = zeros(1,n);
        y(i) = 1;
        y1ddash = 0;
        ynddash = 0;
        interpol_matrix(1:4*(n-1),i) = cubic_sqrt_spline_recipe(x,y,y1ddash,ynddash);
    end
end

%the coefficient matrix is the same for g and for h. But there are also
%infinity conditions.
%We power series expand g, h around infinity:
%g = cx+d+e/x+f/x^2, . . . . . h=bx^2+cx+d+e/x+f/x^2

coefficient_matrix = zeros(8*(n-1)+9,2*n+4);
%the first n+2 columns are dedicated to g
%second n+2 columns dedicated to h
coefficient_matrix(1:4*(n-1),1:n+2) = interpol_matrix;
coefficient_matrix(4*(n-1)+1:8*(n-1),n+3:2*n+4) = interpol_matrix;

%fits cx+d+e/x+f/x^2 through gn,gn',gn'',gn-1
g_infty_coeffs = zeros(4,n+2);
g_infty_coeffs(1,n) = 1;
g_infty_coeffs(2,1:n+2) = [3*x(n)^2, 2*x(n), 1, 0]*interpol_matrix(4*(n-2)+1:4*(n-2)+4,:);
g_infty_coeffs(3,n+2) = 1;
g_infty_coeffs(4,n-1) = 1;
g_infty_coeffs = [x(n), 1, x(n)^(-1), x(n)^(-2); 1, 0, -x(n)^(-2), -2*x(n)^(-3); ...
    0, 0, 2*x(n)^(-3), 6*x(n)^(-4); x(n-1), 1, x(n-1)^(-1), x(n-1)^(-2)]^(-1)*g_infty_coeffs;

%g_infty_coeffs = zeros(3,n+2);
%g_infty_coeffs(1,n) = 1;
%g_infty_coeffs(2,n-1) = 1;
%g_infty_coeffs(3,n-2) = 1;
%g_infty_coeffs = [x(n),1,x(n)^(-1); x(n-1),1,x(n-1)^(-1); ...
%    x(n-2),1,x(n-2)^(-1)]\g_infty_coeffs;

coefficient_matrix(8*(n-1)+1:8*(n-1)+4,1:n+2) = g_infty_coeffs;


%fits bx^2+cx+d+e/x+f/x^2 through hn,hn',hn'',hn-1,hn-2
h_infty_coeffs = zeros(5,n+2);
%h_infty_coeffs(1,n) = 1;
%h_infty_coeffs(2,n-1) = 1;
%h_infty_coeffs(3,n-2) = 1;
%h_infty_coeffs(4,n-3) = 1;
%h_infty_coeffs = [x(n)^2,x(n),1,x(n)^(-1); x(n-1)^2,x(n-1),1,x(n-1)^(-1); ...
%    x(n-2)^2,x(n-2),1,x(n-2)^(-1); x(n-3)^2,x(n-3),1,x(n-3)^(-1)]\h_infty_coeffs;

h_infty_coeffs(1,n) = 1;
h_infty_coeffs(2,1:n+2) = [3*x(n)^2, 2*x(n), 1, 0]*interpol_matrix(4*(n-2)+1:4*(n-2)+4,:);
h_infty_coeffs(3,n+2) = 1;
h_infty_coeffs(4,n-1) = 1;
h_infty_coeffs(5,n-2) = 1;
h_infty_coeffs = [x(n)^2, x(n), 1, x(n)^(-1), x(n)^(-2); 2*x(n), 1, 0, -x(n)^(-2), -2*x(n)^(-3); ...
    2, 0, 0, 2*x(n)^(-3), 6*x(n)^(-4); x(n-1)^2, x(n-1), 1, x(n-1)^(-1), x(n-1)^(-2); ...
    x(n-2)^2, x(n-2), 1, x(n-2)^(-1), x(n-2)^(-2)]^(-1)*h_infty_coeffs;
coefficient_matrix(8*(n-1)+5:8*(n-1)+9,n+3:2*n+4) = h_infty_coeffs;

%creates the matrix!!
%creates the matrix to measure all the pressures
kernel_matrix = zeros(2*n,8*(n-1)+9);
for j = 1:n
    [pressurerow, shearrow] = pressure_cubic_recipe(z(j),x);
    kernel_matrix(j, :) = pressurerow;
    kernel_matrix(n+j,:) = shearrow;
end

%sets conditions at infinity: in order
%x coeff of g, 1 coeff of g
%x^2 coeff of h, x coeff of h
g_infty_conditions = g_infty_coeffs(1:2,:);
h_infty_conditions = h_infty_coeffs(1:2,:);

equation_matrix = zeros(2*n+4,2*n+4);
equation_matrix(1:2*n,1:2*n+4) = kernel_matrix*coefficient_matrix;
equation_matrix(2*n+1:2*n+2,1:n+2) = g_infty_conditions;
equation_matrix(2*n+3:2*n+4,n+3:2*n+4) = h_infty_conditions;

inhomog = zeros(2*n+4,1);
inhomog(2*n+1) = 0;              %x coeff of g
inhomog(2*n+2) = (-Pload+6*Mload);       %1 coeff of g
inhomog(2*n+3) = 0;               %x^2 coeff of h
inhomog(2*n+4) = (12*Mload);        %x coeff of h

%normalizes it / dimensionalizes
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

%finds the solution!!!
soln_vector = equation_matrix\inhomog;
conditioning = rcond(equation_matrix);

%finds KI, KII
KI = interpol_matrix(4,1:n+2)*soln_vector(n+3:2*n+4);
KII = interpol_matrix(4,1:n+2)*soln_vector(1:n+2);

measure = tan((0.05:0.1:n-0.05)*atan(endpoint)/n);
m = length(measure);
pressure = zeros(1,m);
shear = zeros(1,m);

for i = 1:m
    [pressurerow,shearrow] = pressure_cubic_recipe(measure(i),x);
    pressure(i) = pressurerow*coefficient_matrix*soln_vector;
    shear(i) = shearrow*coefficient_matrix*soln_vector;
end
        
        
        
        
        