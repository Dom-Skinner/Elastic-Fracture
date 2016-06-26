%uses a proper cubic spline (as suggested by Lister not Gunnar)
%does NOT use an interpol matrix (omg i was so stupid)
%uses representation x^(-1/2) * (ax^3 + bx^2 + cx + d) up to panel t
%uses representation (bx^2 + cx + d + e/x) past panel t, to infinity
%uses representation from first panel to crack tip
%EXCEPT for g, uses a representation cx + d + e/x + f/x^2 ONLY IN THE LAST
%PANEL WHICH GOES TO INFINITY

%total of 8*n unknowns
%measure pressure at k points: 2*k equations
%specify continuity at n-1 points: 6*(n-1) equations
%additional 4 equations for far field.
%so need k = n+1 points to measure.

n = 50;
endpoint = 20;
Pload = 0;
Mload = 1;

%vector of panel joining points: we will only use n-1 of them though
x = tan((1:n)*(atan(endpoint))/(n+1));

%vector of measurement points
z = tan((0.5:1:n+0.5)*(atan(endpoint))/(n+1));

t = round(n/2);

%we want our solution to be of the form: a vector of 8*(n-1) coefficients.
%two batches of 4*(n-1) coefficients, in turn come in batches of 4 for
%ax^3+bx^2+cx+d.

%creates the matrix to measure all the pressures
kernel_matrix = zeros(2*(n+1),8*n);
for j = 1:n+1
    [pressurerow, shearrow] = pressure_cubic(z(j),x);
    kernel_matrix(j, :) = pressurerow;
    kernel_matrix(n+1+j,:) = shearrow;
end



%creates a matrix to specify continuity
cont_matrix = zeros(n-1, 4*n);
%for first half of panels
for i=1:t-1
    cont_matrix(i, 4*(i-1)+1:4*i) = -[x(i)^3,x(i)^2,x(i),1];
    cont_matrix(i, 4*i+1:4*(i+1)) = [x(i)^3,x(i)^2,x(i),1];
end
%for the joining point
cont_matrix(t, 4*(t-1)+1:4*t) = -[x(t)^(5/2), x(t)^(3/2), x(t)^(1/2), x(t)^(-1/2)];
cont_matrix(t, 4*t+1:4*(t+1)) = [x(t)^2, x(t), 1, x(t)^(-1)];
%for the second half of the panels
for i=t+1:n-1
    cont_matrix(i, 4*(i-1)+1:4*i) = -[x(i)^2,x(i)^1,1,x(i)^(-1)];
    cont_matrix(i, 4*i+1:4*(i+1)) = [x(i)^2,x(i)^1,1,x(i)^(-1)];
end



%creates a matrix to specify the continuity of first derivative
first_deriv_cont = zeros(n-1, 4*n);
%for first half of panels
for i = 1:t-1
    first_deriv_cont(i, 4*(i-1)+1:4*i) = -[3*x(i)^2, 2*x(i), 1, 0];
    first_deriv_cont(i, 4*i+1:4*(i+1)) = [3*x(i)^2, 2*x(i), 1, 0];
end
%for the joining panel
first_deriv_cont(t, 4*(t-1)+1:4*t) = -[(5/2)*x(t)^(3/2), (3/2)*x(t)^(1/2), ...
    (1/2)*x(t)^(-1/2), -(1/2)*x(t)^(-3/2)];
first_deriv_cont(t, 4*t+1:4*(t+1)) = [2*x(t), 1, 0, -x(t)^(-2)];
%for the remaining panels
for i = t+1:n-1
    first_deriv_cont(i, 4*(i-1)+1:4*i) = -[2*x(i), 1, 0, -x(i)^(-2)];
    first_deriv_cont(i, 4*i+1:4*(i+1)) = [2*x(i), 1, 0, -x(i)^(-2)];
end


%creates a matrix to specify continuity of second derivative
second_deriv_cont = zeros(n-1, 4*n);
%for first half of panels
for i = 1:t-1
    second_deriv_cont(i, 4*(i-1)+1:4*i) = -[6*x(i), 2, 0, 0];
    second_deriv_cont(i, 4*i+1:4*(i+1)) = [6*x(i), 2, 0, 0];
end
%for the joining panel
second_deriv_cont(t, 4*(t-1)+1:4*t) = -[(5/2)*(3/2)*x(t)^(1/2), (3/2)*(1/2)*x(t)^(-1/2), ...
    (1/2)*(-1/2)*x(t)^(-3/2), (1/2)*(3/2)*x(t)^(-5/2)];
second_deriv_cont(t, 4*t+1:4*(t+1)) = [2, 0, 0, 2*x(t)^(-3)];
%for the remaining panels
for i = t+1:n-1
    second_deriv_cont(i, 4*(i-1)+1:4*i) = -[2, 0, 0, 2*x(i)^(-3)];
    second_deriv_cont(i, 4*i+1:4*(i+1)) = [2, 0, 0, 2*x(i)^(-3)];
end  

g_cont_matrix = cont_matrix;
h_cont_matrix = cont_matrix;

g_first_deriv_cont = first_deriv_cont;
h_first_deriv_cont = first_deriv_cont;

g_second_deriv_cont = second_deriv_cont;
h_second_deriv_cont = second_deriv_cont;

%modifies the g condition to fit the representation at infinity
%function itself
g_cont_matrix(n-1,4*(n-2)+1:4*(n-1)) = -[x(n-1)^2, x(n-1), 1, x(n-1)^(-1)];
g_cont_matrix(n-1,4*(n-1)+1:4*n) = [x(n-1),1,x(n-1)^(-1),x(n-1)^(-2)];
%first derivative
g_first_deriv_cont(n-1, 4*(n-2)+1:4*(n-1)) = -[2*x(n-1), 1, 0, -x(n-1)^(-2)];
g_first_deriv_cont(n-1, 4*(n-1)+1:4*n) = [1, 0, -x(n-1)^(-2), -2*x(n-1)^(-3)];
%second derivative
g_second_deriv_cont(n-1, 4*(n-2)+1:4*(n-1)) = -[2, 0, 0, 2*x(i)^(-3)];
g_second_deriv_cont(n-1, 4*(n-1)+1:4*n) = [0, 0, 2*x(n-1)^(-3), 6*x(n-1)^(-4)];



%creates a matrix of equations.
%first 4(n-1) rows: the pressures and shears
%next 2(n-2) rows: continuity of g
%next 2(n-2) rows: continuity of h
%final 4 rows: far field conditions on g (x then 1 coefficient)
%and on h (x^2 then x coefficient)

equation_matrix = zeros(8*n,8*n);

equation_matrix(1:2*(n+1),:) = kernel_matrix;

equation_matrix(2*(n+1)+1:2*(n+1)+(n-1), 1:4*n) = g_cont_matrix;
equation_matrix(2*(n+1)+(n-1)+1:2*(n+1)+2*(n-1), 1:4*n) = g_first_deriv_cont;
equation_matrix(2*(n+1)+2*(n-1)+1:2*(n+1)+3*(n-1), 1:4*n) = g_second_deriv_cont;

equation_matrix(2*(n+1)+3*(n-1)+1:2*(n+1)+4*(n-1), 4*n+1:8*n) = h_cont_matrix;
equation_matrix(2*(n+1)+4*(n-1)+1:2*(n+1)+5*(n-1), 4*n+1:8*n) = h_first_deriv_cont;
equation_matrix(2*(n+1)+5*(n-1)+1:2*(n+1)+6*(n-1), 4*n+1:8*n) = h_second_deriv_cont;

%for x coeff of g
equation_matrix(8*n-3, 4*(n-1)+1) = 1;
%for 1 coeff of g
equation_matrix(8*n-2, 4*(n-1)+2) = 1;
%for x^2 coeff of h
equation_matrix(8*n-1, 4*n+4*(n-1)+1) = 1;
%for the x coeff of h
equation_matrix(8*n, 4*n+4*(n-1)+2) = 1;


%the inhomogenous forcing
inhomog = zeros(8*n,1);
inhomog(8*n-3) = 0;       %x coeff of g
inhomog(8*n-2) = (-Pload+6*Mload);       %1 coeff of g
inhomog(8*n-1) = 0;       %x^2 coeff of h
inhomog(8*n) = (12*Mload);         %x coeff of h

%normalization/dimensionalization factor
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

soln_coeffs = equation_matrix\inhomog;

KI = soln_coeffs(4*n+4);
KII = soln_coeffs(4);



