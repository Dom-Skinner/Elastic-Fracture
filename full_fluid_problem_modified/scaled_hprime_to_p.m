%this function takes a given h' and returns two the pressure map, and the
%derivative of the pressure map.

function [p, dp] = scaled_hprime_to_p(x,z,coeffs,lambda, h_coefficient_matrix)

n = length(x);

%input is a 4*n vector of g' and h'. coeffs

%finds the relevant coefficients of h

%h_coefficient_matrix = hprime_to_h(x);
h_coeffs = h_coefficient_matrix*coeffs(2*n+1:4*n);

%finds the pressure: an n-1 vector
pressure = pprime_to_p(x,z,h_coeffs);
%finds the pressure differential: a (n-1)xn matrix
deriv_pressure = pressure_map_derivative(x,z,h_coeffs,h_coefficient_matrix);

p = zeros(4*n,1);
dp = zeros(4*n,4*n);

p(1:n-1) = lambda*pressure;
dp(1:n-1,2*n+1:4*n) = lambda*deriv_pressure;

%also sets the M and P BC's inside p

adjust_pen = 1-2*lambda/(3*x(n-1)) - 4*lambda^2 *log(x(n-1)) / x(n-1)^2 ;
adjust_end = 1-2*lambda/(3*x(n)) - 4*lambda^2 *log(x(n)) / x(n)^2 ;
p(4*n-3) = 0.5*adjust_pen;
p(4*n-2) = 0.5*adjust_end;
p(4*n-1) = adjust_pen;
p(4*n)   = adjust_end;

return
end