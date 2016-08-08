function [theta,h] = linear_perturbation_solve(n,t,x,z, h0,h0_z,l0,kernel_matrix,s)
% Solves the linear perturbation problem.


%the data points

% We will have kernel_matrix as an input, since it takes so long to compute
%finds the appropriate elasticity kernel
%[kernel_matrix, ~] = pressure_shear_matrix_s(x,z,s);


%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h_s(x,s,t);
%R = lubrication_integral_s(x,z,n,t,h_coefficient_matrix,h0,l0,0.5);
R = zeros(2*(n-1),2*n);
R(1:n-1,n+1:2*n) = lubrication_integral(x,z,n,t,h_coefficient_matrix,h0,h0_z,l0,s);

%BT = kernel_matrix;
rescale = diag([z.^(1-s), ones(1,n-1)]); % To account for singularity in p
A = zeros(2*n,2*n);
A(1:2*(n-1),:) = rescale*(kernel_matrix-R);
%the limit of g' at infinity
A(2*n-1,n) = 1;
%the limit of h'' at infinity
A(2*n,2*n-1) = -1/(x(n)-x(n-1));
A(2*n,2*n) = 1/(x(n)-x(n-1));

c = zeros(2*n,1);
c(2*n-1,1) = 0.5;
c(2*n,1) = 1;

theta = A\c;

h_coeffs = h_coefficient_matrix*theta(n+1:2*n,1);
h = h_eval(h_coeffs,x,n,t,s);
end