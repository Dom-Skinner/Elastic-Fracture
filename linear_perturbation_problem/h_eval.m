function h = h_eval(h_coeffs,x,n,t,s)
% Simply evaluates h, given its coefficients
h = zeros(n,1);

h(1: t-1) = x(1:t-1)'.^(1+s) .* h_coeffs(1:3:3*(t-1)-2) + x(1:t-1)'.^(s) ...
    .*h_coeffs(2:3:3*(t-1)-1) + h_coeffs(3:3:3*(t-1));
h(t:n) = x(t:n)'.^2 .* h_coeffs(3*t-2:3:3*n) + x(t:n)' ...
    .*h_coeffs(3*t-1:3:3*n) + h_coeffs(3*t:3:3*n);


end