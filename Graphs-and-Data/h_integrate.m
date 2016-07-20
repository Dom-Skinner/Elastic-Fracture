function h = h_integrate(h_prime,x,n,t,h_coefficient_matrix)

h_coeffs = h_coefficient_matrix*h_prime;

h = zeros(n,1);

for k = 1:t-1
    h(k) = x(k).^(1.5)*h_coeffs(3*k-2) + x(k)^(0.5)*h_coeffs(3*k-1)+...
        h_coeffs(3*k);
end

for k = t:n
    h(k) = x(k)^2*h_coeffs(3*k-2) + x(k)*h_coeffs(3*k-1)+...
        h_coeffs(3*k);
end


end