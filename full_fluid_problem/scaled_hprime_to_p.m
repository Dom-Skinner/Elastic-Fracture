%this function takes a given h' and returns two the pressure map, and the
%derivative of the pressure map.

function [p, dp] = scaled_hprime_to_p(x,z,hprime,lambda, h_coefficient_matrix,t)

n = length(x);

%input is a 2*n vector of g' and h'.

%finds the relevant coefficients of h

%h_coefficient_matrix = hprime_to_h(x);
h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);

%finds the pressure: an n-1 vector
pressure = pprime_to_p(x,z,h_coeffs,t);
%finds the pressure differential: a (n-1)xn matrix
deriv_pressure = pressure_map_derivative(x,z,h_coeffs,h_coefficient_matrix,t);

%finds the adjustments to things
%[adjust, d_adjust] = bending_p_adjust(x,h_coeffs,h_coefficient_matrix);

p = zeros(2*n,1);
dp = zeros(2*n,2*n);

p(1:n-1) = lambda*pressure;
dp(1:n-1,n+1:2*n) = lambda*deriv_pressure;

%also sets the M and P BC's inside p

%the P condition


adjust= 1-2*lambda./(3*x) - 4*lambda^2 *log(x) ./ x.^2 ;

p(2*n-1) = 0.5*adjust(n);
p(2*n) = adjust(n-1);
return
end