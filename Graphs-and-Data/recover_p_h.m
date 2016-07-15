% function that recovers p and h given hprime.
% Probably not the most efficient, but it works
% Notably it will give the lambda scaled out versions of h and p
% AKA h-tilde and p-tilde

function [p, h] = recover_p_h(lambda,x,z,num,hprime_data,n,t)
lend = numel(lambda);

h = zeros(numel(num),lend);
p = zeros(numel(num),lend);
h_coefficient_matrix = hprime_to_h(x);
[kernel_matrix, ~] = pressure_shear_matrix(x,z);


for l = 1:lend
    hprime = hprime_data(:,l);
    h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);
    h_temp = zeros(1,t);
    for r = 1:t-1
        h_temp(r) = h_coeffs(3*r -2)*x(r)^(3/2) + ...
            h_coeffs(3*r-1)*x(r)^(1/2) + h_coeffs(3*r);
    end
    h(:,l) = lambda(l).^(-1/3).*h_temp(num);
    pressure = kernel_matrix*hprime;
    p(:,l) = lambda(l).^(-1/3).*pressure(num);
    
end

end