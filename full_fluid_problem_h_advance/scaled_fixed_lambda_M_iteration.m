function [KI,KII, hprime_new,p] = scaled_fixed_lambda_M_iteration(x,z,v,w,t,xmax, lambda, tol, hprime_start)
% Function that iterates, finding K for some given value of lambda,
% up to some tolerance, tol. 

% Some typical values previously used
%lambda = 0.02, n = 100, tol = 10^(-8), xmax = 10;


% For now ... 

nx = numel(x);
nv = numel(v);
vt = t + nv-nx;
L = v(1+nv-nx);

%finds the appropriate elasticity kernel
%K11 = pressure_matrix_g(v,w(end-nx+2:end),vt,lambda);
%K12 = pressure_matrix_h(x,z,t,lambda);
K11 = pressure_matrix_g(v,w,vt,lambda);
K12 = pressure_matrix_h(x,w-L,t,lambda);
K21 = shear_matrix_g(v,w,vt,lambda);
K22 = shear_matrix_h(x,w-L,t,lambda);

kernel_matrix = [K11(1+nv-nx:end,:),K12(1+nv-nx:end,:) ; K21,K22];


%in short, we wish to solve Ah = f(h), where f is our map from h' to p.

%to do this, we try find dh such that
%A(h+dh) = f(h+dh), or Ah + Adh = f(h) + df(h)*dh
%or dh = (A-df(h))^(-1) * (f(h)-Ah), a version of newton iteration
%a.k.a. if h is the first iterate,
%then h + (A-df(h))^(-1) * (f(h)-Ah) is the next iterate.
%matrix A is kernel_matrix as given above, however we need to add some
%boundary conditions:
%these will be inside A and f, in the form of extra equations.

A = zeros(nx+nv,nx+nv);
A(1:nx+nv-2,:) = kernel_matrix;

%the limit of g' at infinity
A(end-1,nv) = 1;

%the limit of h'' at infinity
A(end,end-1:end) = [-1/(x(end)-x(end-1)),1/(x(end)-x(end-1))];

conditioning = rcond(A);

%saves a matrix to convert h' to h
h_coefficient_matrix = hprime_to_h(x,t,lambda);

fprintf('\n n = %d xmax = %d lambda = %6.4g L = %4.4g\n', nv, xmax, lambda,L)
fprintf(' condition number = %6.4e \n', conditioning)
fprintf('iteration,  KI,       dKI,         KII,          dKII           global diff\n');

hprime_new = hprime_start;
KI_new = hprime_new(nv+1);
KII_new = hprime_new(1);
iteration_number = 0;
KI_diff = tol+1; % Just to ensure that it iterates at least once

while abs(KI_diff) >= tol && iteration_number<25
    hprime_old = hprime_new;
    KI_old = KI_new;
    KII_old = KII_new;

    [p,dp] = scaled_hprime_to_p(x,z,hprime_old(nv+1:end),lambda,h_coefficient_matrix,nv+nx,t);
    hprime_new = hprime_old + (A-dp)\(p-A*hprime_old);
    %p1 = K11*hprime_new(1:nv);
    %p2 = K12*hprime_new(nv+1:end);
    %p = [K11, K12]*hprime_new;
    
    KI_new = hprime_new(nv+1);
    KII_new = hprime_new(1);
    KI_diff = KI_new - KI_old;
    KII_diff = KII_new - KII_old;
    global_diff = norm(hprime_new - hprime_old);
    
    iteration_number = iteration_number + 1;
    fprintf('%3d & %11.4e & %11.4e & %11.4e & %11.4e & %11.4e \\\\ \n',iteration_number,...
        KI_new, KI_diff,KII_new, KII_diff, global_diff);
   
end

KI  = KI_new;
KII = KII_new;
if iteration_number == 25
    fprintf('Failure to coverge\n');
end

return
end