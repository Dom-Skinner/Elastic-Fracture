%this function takes a given h' and returns two the pressure map, and the
%derivative of the pressure map.

function [p, dp] = hprime_to_p_num(x,z,hprime,lambda, Mload, Pload, h_coefficient_matrix)

n = length(x);

%input is a 2*n vector of g' and h'.

%finds the relevant coefficients of h

%h_coefficient_matrix = hprime_to_h(x);
h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);

%finds the pressure: an n-1 vector
pressure = pprime_to_p(x,z,h_coeffs,lambda);
%finds the pressure differential: a (n-1)xn matrix
deriv_pressure_num = pressure_map_derivative_num(x,z,h_coeffs,lambda,h_coefficient_matrix);

p = zeros(2*n,1);
dp = zeros(2*n,2*n);

p(1:n-1) = pressure;
dp(1:n-1,n+1:2*n) = deriv_pressure_num;

%also sets the M and P BC's inside p

%the P condition
p(2*n-1) = (-Pload+6*Mload)*(2*pi)^(3/2)/(8*pi);
p(2*n) = (12*Mload)*(2*pi)^(3/2)/(8*pi);

return
end