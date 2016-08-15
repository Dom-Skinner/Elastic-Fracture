function gprime = h_to_g(x,z,v,w,t,hprime,lambda)
% Converts hprime into gprime.
% assumes hprime is literally just hprime and not [g',h'] as it sometimes
% is written.

nx = numel(x);
nv = numel(v);
vt = t + nv-nx;
L = v(1+nv-nx);

%finds the appropriate elasticity kernel
%K11 = pressure_matrix_g(v,w(end-nx+2:end),vt,lambda);
%K12 = pressure_matrix_h(x,z,t,lambda);

K21 = shear_matrix_g(v,w,vt,lambda);
K22 = shear_matrix_h(x,w-L,t,lambda);

A = zeros(nv,nv);
B = zeros(nv,nx);
A(1:end-1,:) =  K21;
B(1:end-1,:) = K22;

%the limit of g' at infinity
A(end,end) = 1;
%the limit of h'' at infinity

c =  - B*hprime;
c(end) = 0.5*( 1-2*lambda./(3*x(end)) - 4*lambda^2 *log(x(end)) ./ x(end).^2 ) ;

gprime = A\c;
end