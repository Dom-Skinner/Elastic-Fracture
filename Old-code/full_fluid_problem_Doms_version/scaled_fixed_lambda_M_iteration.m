function [KI, hprime_new,p] = scaled_fixed_lambda_M_iteration(x,z,n,t,xmax, lambda, tol, hprime_start)
% Function that iterates, finding K for some given value of lambda,
% up to some tolerance, tol. ~Dom

% Some typical values previously used
%lambda = 0.02, n = 100, tol = 10^(-8), xmax = 10;

%the data points


%finds the appropriate elasticity kernel
[kernel_matrix, interpolate_matrix] = pressure_shear_matrix(x,z,t,lambda);

%in short, we wish to solve Ah = f(h), where f is our map from h' to p.
%to do this, we try find dh such that
%A(h+dh) = f(h+dh), or Ah + Adh = f(h) + df(h)*dh
%or dh = (A-df(h))^(-1) * (f(h)-Ah), a version of newton iteration
%a.k.a. if h is the first iterate,
%then h + (A-df(h))^(-1) * (f(h)-Ah) is the next iterate.
%matrix A is kernel_matrix as given above, however we need to add some
%boundary conditions:
%these will be inside A and f, in the form of extra equations.

A = zeros(2*n,2*n);
A(1:2*(n-1),:) = kernel_matrix;
%the limit of g' at infinity

%the limit of h'' at infinity

A(2*n-1,n) = 1;
A(2*n,:) = interpolate_matrix(3*n-1,:);

conditioning = rcond(A);

%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h(x,t,lambda);

fprintf('\n n = %d xmax = %d lambda = %6.4g \n', n, xmax, lambda)
fprintf(' condition number = %6.4e \n', conditioning)
fprintf('iteration,  K,       dK,          global diff\n');

hprime_new = hprime_start;
K_new = hprime_new(n+1);
iteration_number = 0;
K_diff = tol+1; % Just to ensure that it iterates at least once

while abs(K_diff) >= tol && iteration_number<13
    hprime_old = hprime_new;
    K_old = K_new;
    

    [p,dp] = scaled_hprime_to_p(x,z,hprime_old,lambda,h_coefficient_matrix,t);
    hprime_new = hprime_old + (A-dp)\(p-A*hprime_old);
    
    K_new = hprime_new(n+1);
    K_diff = K_new - K_old;
    global_diff = norm(hprime_new - hprime_old);
    
    iteration_number = iteration_number + 1;
    fprintf('%3d & %11.4e & %11.4e & %11.4e \\\\ \n',iteration_number, K_new, K_diff, global_diff);
   
end

KI = K_new;
if iteration_number == 13
    fprintf('Failure to coverge');
    KI=-1;
end

return
end