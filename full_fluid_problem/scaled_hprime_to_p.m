%this function takes a given h' and returns two the pressure map, and the
%derivative of the pressure map.

function [p, dp] = scaled_hprime_to_p(x,z,hprime,lambda, h_coefficient_matrix)

n = length(x);

%input is a 2*n vector of g' and h'.

%finds the relevant coefficients of h

%h_coefficient_matrix = hprime_to_h(x);
h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);

%finds the pressure: an n-1 vector
pressure = pprime_to_p(x,z,h_coeffs);
%finds the pressure differential: a (n-1)xn matrix
deriv_pressure = pressure_map_derivative(x,z,h_coeffs,h_coefficient_matrix);

%finds the adjustments to things
[adjust, d_adjust] = bending_p_adjust(x,h_coeffs,h_coefficient_matrix);

p = zeros(2*n,1);
dp = zeros(2*n,2*n);

p(1:n-1) = lambda*pressure;
dp(1:n-1,n+1:2*n) = lambda*deriv_pressure;

%also sets the M and P BC's inside p

%scale = 12*(2*pi)^(3/2)/(8*pi)*Mload;

%h_doubledash = scale - 2*lambda/(3*scale^2*x(n));

%the P condition
%p(2*n-1) = -Pload*(2*pi)^(3/2)/(8*pi) + 0.5*h_doubledash;
%p(2*n) = h_doubledash;

%the P condition

p(2*n-1) = 0.5 + 0.5*lambda*adjust;
p(2*n) = 1 + lambda*adjust;

dp(2*n-1,n+1:2*n) = 0.5*lambda*d_adjust;
dp(2*n,n+1:2*n) = lambda*d_adjust;




return
end