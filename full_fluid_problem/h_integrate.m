function h = h_integrate(h_prime,x,n,t,h_coefficient_matrix,s)
% Recovers h from h_prime, given the coefficient matrix
h_coeffs = h_coefficient_matrix*h_prime;
h = h_eval(h_coeffs,x,n,t,s);

end