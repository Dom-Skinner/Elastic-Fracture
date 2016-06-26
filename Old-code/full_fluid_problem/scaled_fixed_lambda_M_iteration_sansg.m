function [KI, hprime_new,p] = scaled_fixed_lambda_M_iteration_sansg(n,xmax, lambda, tol, hprime_start)

%lambda = 20;
%tol = 10^(-8);

%number of data points
%n = 100;
t = round(n/2);
%xmax = 10;

%the data points
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
%x = zeros(1,n);
%x(1:t) = tan((0:t-1)*(pi/4)/(t-1)).^2;
%x(t+1:n) = (1:n-t)*(xmax-1)/(n-t) + ones(1,n-t);

%x = zeros(1,n);
%x(1:t) = tan((0:t-1)*(pi/4)/(t-1)).^2;
%x(t+1:n) = tan((1:n-t)*(atan(xmax)-atan(1))/(n-t));

z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
%z = zeros(1,n-1);
%z(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1)).^2;
%z(t:n-1) = (0.5:n-t-0.5)*(xmax-1)/(n-t) + ones(1,n-t);

%z = zeros(1,n-1);
%z(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1)).^2;
%z(t:n-1) = tan((0.5:n-t-0.5)*(atan(xmax)-atan(1))/(n-t));

%hprime_start = zeros(2*n,1);
%hprime_start(1:n) = ones(n,1);
%hprime_start(n+1:2*n) = x' + 1;

%finds the appropriate elasticity kernel
[kernel_matrix, interpolate_matrix] = pressure_shear_matrix_sansg(x,z);

%in short, we wish to solve Ah = f(h), where f is our map from h' to p.
%to do this, we try find dh such that
%A(h+dh) = f(h+dh), or Ah + Adh = f(h) + df(h)*dh
%or dh = (A-df(h))^(-1) * (f(h)-Ah), a version of newton iteration
%a.k.a. if h is the first iterate,
%then h + (A-df(h))^(-1) * (f(h)-Ah) is the next iterate.
%matrix A is kernel_matrix as given above, however we need to add some
%boundary conditions:
%these will be inside A and f, in the form of extra equations.

A = zeros(n,n);
A(1:n-1,:) = kernel_matrix;
%the limit of h'' at infinity
A(n,:) = interpolate_matrix(n,:);

conditioning = rcond(A);

%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h_sansg(x);

hprime_old = hprime_start;
K_old = hprime_old(1);

[p,dp] = scaled_hprime_to_p_sansg(x,z,hprime_old,lambda,h_coefficient_matrix);
hprime_new = hprime_old + (A-dp)\(p-A*hprime_old);

K_new = hprime_new(1)*3*sqrt(2*pi);
K_diff = K_new-K_old;
global_diff = norm(hprime_new - hprime_old);

iteration_number = 1;
fprintf('\n n=%7d xmax=%7d lambda=%7d \n', n, xmax, lambda)
fprintf('condition number = %13.7e \n', conditioning)
fprintf('%7d & %13.7e & %13.7e & %13.7e \\\\ \n',iteration_number, K_new, K_diff, global_diff);

while abs(K_diff) >= tol
    hprime_old = hprime_new;
    K_old = K_new;
    
    [p,dp] = scaled_hprime_to_p_sansg(x,z,hprime_old,lambda,h_coefficient_matrix);
    hprime_new = hprime_old + (A-dp)\(p-A*hprime_old);
    
    K_new = hprime_new(1)*3*sqrt(2*pi);
    K_diff = K_new - K_old;
    global_diff = norm(hprime_new - hprime_old);
    
    iteration_number = iteration_number + 1;
    fprintf('%7d & %13.7e & %13.7e & %13.7e \\\\ \n',iteration_number, K_new, K_diff, global_diff);

end

KI = K_new;

return
end