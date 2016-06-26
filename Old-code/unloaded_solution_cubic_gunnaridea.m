%function [KI, KII, soln, x, equation_matrix, kernel_matrix, interpolate_matrix] = ...
%unloaded_solution_cubic(n, endpoint, Pload,Mload)

n = 400;
endpoint = 50;
Pload = 0;
Mload = 1;
    
%we use x^(-1/2)*(ax^3+bx^2+cx+d) spline for first t panels, inc. tip
%use (cx+d+e/x+f/x^2) spline for the final panels, inc. infinity

%this means has to measure pressure at twice number of points as otherwise

%does the spline up to endpoint
%distributes points using tan

%vector of h' points
x = tan((1:n)*(atan(endpoint))/n);

%vector of pressure points
outpoint = tan((1.25:0.5:n-0.25)*(atan(endpoint))/n);

%matrix for the interpolation, has to sort out infinity though: needs to
%add another 9 values, for g: cx+d+e/x+f/x^2, for h: bx^2+...+f/x^2

interpolate_matrix = zeros(8*(n-1)+8,4*n);
interpolate_matrix(1:4*(n-1),1:2*n) = interpolate_cubic(x);
interpolate_matrix(4*(n-1)+1:8*(n-1), 2*n+1:4*n) = interpolate_cubic(x);

g_infty_coefficients = [x(n-1), 1, x(n-1)^(-1), x(n-1)^(-2); ... %g(x(n-1))
    x(n), 1, x(n)^(-1), x(n)^(-2); ... %g(x(n))
    1, 0, -x(n-1)^(-2), -2*x(n-1)^(-3); ... %g'(x(n-1))
    1, 0, -x(n)^(-2), -2*x(n)^(-3)]^(-1); %g'(x(n))

h_infty_coefficients = [x(n-1)^2, x(n-1), 1, x(n-1)^(-1); ... %h(x(n-1))
    x(n)^2, x(n), 1, x(n)^(-1); ... %h(x(n))
    2*x(n-1), 1, 0, -x(n-1)^(-2); ... %h'(x(n-1))
    2*x(n), 1, 0, -x(n)^(-2)]^(-1); %h'(x(n))

interpolate_matrix(8*(n-1)+1:8*(n-1)+4, [n-1,n,2*n-1,2*n]) = g_infty_coefficients;
interpolate_matrix(8*(n-1)+5:8*(n-1)+8, [3*n-1,3*n,4*n-1,4*n]) = h_infty_coefficients;

%creates matrix for the kernels
%has 2 x (2n-2) rows, one set of 2n-2 for pressure, other set for the shear
kernel_matrix = zeros(2*(2*n-2), 8*(n-1)+8);

%fills the rows of Kmatrix
for j = 1:2*n-2
    z = outpoint(j);
    [pressurerow, shearrow] = pressure_cubic_gunnar(z,x);
    kernel_matrix(j, :) = pressurerow;
    kernel_matrix(2*n-2+j,:) = shearrow;
end

pressure_matrix = kernel_matrix*interpolate_matrix;

%needs to form the final 4 equations for the matrix: from boundary
%conditions at infinity.

g_infty_conditions = zeros(2,4*n);
h_infty_conditions = zeros(2,4*n);

g_infty_conditions(:,[n-1, n, 2*n-1, 2*n]) = g_infty_coefficients(1:2,:);
h_infty_conditions(:,[3*n-1, 3*n, 4*n-1, 4*n]) = h_infty_coefficients(1:2,:);


%forms the matrix to be inverted

equation_matrix = zeros(4*n,4*n);
equation_matrix(1:4*n-4,:) = pressure_matrix;
equation_matrix(4*n-3:4*n-2,:) = g_infty_conditions;
equation_matrix(4*n-1:4*n,:) = h_infty_conditions;

%gets the inhomogenous forcing
inhomog = zeros(4*n,1);

inhomog(4*n-3) = 0; %(-6*Nload); %x coeff of g'
inhomog(4*n-2) = (-Pload+6*Mload); %constant coeff of g'
inhomog(4*n-1) = 0; %(6*Nload); %x^2 coeff of h'
inhomog(4*n) = (12*Mload); %x coeff of h'

%adds normalization factor
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

%soln
soln = equation_matrix\inhomog;

KI = (x(2)*soln(2*n+1)-soln(2*n+2)*x(1))/(x(2)-x(1));
KII = (x(2)*soln(1)-soln(2)*x(1))/(x(2)-x(1));
    
    




