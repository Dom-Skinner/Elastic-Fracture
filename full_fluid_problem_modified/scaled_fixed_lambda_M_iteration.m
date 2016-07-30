function [KI, coeff_new,p] = scaled_fixed_lambda_M_iteration(n,xmax, lambda, tol, coeff_start)
% Function that iterates, finding K for some given value of lambda,
% up to some tolerance, tol. 

%the data points
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

%finds the appropriate elasticity kernel
kernel_matrix = pressure_shear_matrix(x,z);

%in short, we wish to solve Ah = f(h), where f is our map from h' to p.
%to do this, we try find dh such that
%A(h+dh) = f(h+dh), or Ah + Adh = f(h) + df(h)*dh
%or dh = (A-df(h))^(-1) * (f(h)-Ah), a version of newton iteration
%a.k.a. if h is the first iterate,
%then h + (A-df(h))^(-1) * (f(h)-Ah) is the next iterate.
%matrix A is kernel_matrix as given above, however we need to add some
%boundary conditions:
%these will be inside A and f, in the form of extra equations.

A = zeros(4*n,4*n);
A(1:2*(n-1),:) = kernel_matrix;
% continuity conditions
A(2*n-1:4*(n-1),:) = continuity_conditions(x);
%the limit of g' at infinity
A(4*n-3,  n-1) = x(n-1);
A(4*n-3,2*n-1) = 1;
A(4*n-2,  n  ) = x(n);
A(4*n-2,2*n  ) = 1;
% The limit of h'' at infinity
A(4*n-1,3*n-1) = 1;
A(4*n  ,3*n  ) = 1;
%
%
conditioning = rcond(A);

%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h(x);

fprintf('\n n = %d xmax = %d lambda = %6.4g \n', n, xmax, lambda)
fprintf(' condition number = %6.4e \n', conditioning)
fprintf('iteration,  K,       dK,          global diff\n');

coeff_new = coeff_start;
K_new = coeff_new(3*n+1);
iteration_number = 0;
K_diff = tol+1; % Just to ensure that it iterates at least once

while abs(K_diff) >= tol && iteration_number<20
    coeff_old = coeff_new;
    K_old = K_new;
    
    [p,dp] = scaled_hprime_to_p(x,z,coeff_old,lambda,h_coefficient_matrix);
    coeff_new = coeff_old + (A-dp)\(p-A*coeff_old);
    
    K_new = coeff_new(3*n+1);
    K_diff = K_new - K_old;
    theta_old = coeff_old(2*n+1:3*n).*x' + coeff_old(3*n+1:end);
    theta_new = coeff_new(2*n+1:3*n).*x' + coeff_new(3*n+1:end);
    %global_diff = norm(coeff_new - coeff_old);
    global_diff = norm(theta_new - theta_old);
    
    iteration_number = iteration_number + 1;
    fprintf('%3d & %11.4e & %11.4e & %11.4e \\\\ \n',iteration_number, K_new, K_diff, global_diff);
   
end

KI = K_new;
if iteration_number == 20
    fprintf('Failure to coverge\n');
    KI=-1;
end

return
end