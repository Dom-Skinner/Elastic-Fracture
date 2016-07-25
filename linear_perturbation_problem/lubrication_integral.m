function R = lubrication_integral(x,z,n,t,h_coeff_matrix_s,h0,h0_z,l0,s)

W_co = zeros(n,1);
E_co = zeros(n,1);
R_co = zeros(n,1);

for k = 2:t-1
    [W_co(k),E_co(k),R_co(k) ] = ...
        get_integral_coeff(x(k),x(k+1),h0(k),h0(k+1),s,2/3,l0);
end

for k = t:n-1
    [W_co(k),E_co(k),R_co(k) ] = ...
        get_integral_coeff(x(k),x(k+1),h0(k),h0(k+1),1,2,l0);
end

Infty = 200;
[W_co(n),E_co(n),R_co(n) ] = ...
        get_integral_coeff(x(n),Infty,h0(n),0.5*Infty^2,1,2,l0);

Q = zeros(n-1,3*n);

for k = 1:n-1
    if k < t
        [Q(k,3*k-2),Q(k,3*k-1), Q(k,3*k)] = ...
            get_integral_coeff(z(k),x(k+1),h0_z(k),h0(k+1),s,2/3,l0);
    else
        [Q(k,3*k-2),Q(k,3*k-1), Q(k,3*k)] = ...
            get_integral_coeff(z(k),x(k+1),h0_z(k),h0(k+1),1,2,l0);
    end
    
    for r = k+1: n
        Q(k,3*r-2) = W_co(r);
        Q(k,3*r-1) = E_co(r);
        Q(k,3*r  ) = R_co(r);
    end
    
end

R = Q*h_coeff_matrix_s;

end