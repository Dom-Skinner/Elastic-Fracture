function [theta,h] = scaled_fixed_lambda_M_iteration(n,t,xmax, h0,l0)
% Solves the linear perturbation problem.
% TODO - rename this function to something more sensible. 

%the data points
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

%finds the appropriate elasticity kernel
[kernel_matrix, ~] = pressure_shear_matrix(x,z);
% so this is the part that is exactly as before
%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h(x);
R = lubrication_integral(x,z,h0,l0,h_coefficient_matrix);

%BT = kernel_matrix;

A = zeros(2*n,2*n);
A(1:2*(n-1),:) = kernel_matrix-R;
%the limit of g' at infinity
A(2*n-1,n) = 1;
%the limit of h'' at infinity
A(2*n,2*n-1) = -1/(x(n)-x(n-1));
A(2*n,2*n) = 1/(x(n)-x(n-1));

conditioning = rcond(A);

fprintf(' condition number = %6.4e \n', conditioning)

c = zeros(2*n,1);
c(2*n-1,1) = 0.5;
c(2*n,1) = 1;

theta = A\c;

h_coeffs = h_coefficient_matrix*theta(n+1:2*n,1);
h(1: t-1) = x.^(1.5) .* h_coeffs(1:3:3*(t-1)-2) + x.^(0.5) ...
    .*h_coeffs(2:3:3*(t-1)-1) + h_coeffs(3:3:3*(t-1));
h(t:n) = x.^2 .* h_coeffs(1:3:3*(t-1)-2) + x ...
    .*h_coeffs(2:3:3*(t-1)-1) + h_coeffs(3:3:3*(t-1));
end