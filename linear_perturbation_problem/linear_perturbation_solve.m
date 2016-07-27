function [theta,h,kernel_matrix] = linear_perturbation_solve(n,t,xmax, h0,h0_z,l0)
% Solves the linear perturbation problem.
s = 0.138673;
%s = 0.5;

%the data points
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

%finds the appropriate elasticity kernel
[kernel_matrix, ~] = pressure_shear_matrix_s(x,z,s);
% so this is the part that is exactly as before
%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h_s(x,s);
%R = lubrication_integral_s(x,z,n,t,h_coefficient_matrix,h0,l0,0.5);
R = zeros(2*(n-1),2*n);
R(1:n-1,n+1:2*n) = lubrication_integral(x,z,n,t,h_coefficient_matrix,h0,h0_z,l0,s);

%BT = kernel_matrix;

A = zeros(2*n,2*n);
A(1:2*(n-1),:) = kernel_matrix-R;
%the limit of g' at infinity
A(2*n-1,n) = 1;
%the limit of h'' at infinity
A(2*n,2*n-1) = -1/(x(n)-x(n-1));
A(2*n,2*n) = 1/(x(n)-x(n-1));

%conditioning = rcond(kernel_matrix);

%fprintf(' condition number = %6.4e \n', conditioning)

c = zeros(2*n,1);
c(2*n-1,1) = 0.5;
c(2*n,1) = 1;

theta = A\c;

h_coeffs = h_coefficient_matrix*theta(n+1:2*n,1);
h = h_eval(h_coeffs,x,n,t,s);
end