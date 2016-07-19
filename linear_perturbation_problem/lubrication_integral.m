function R = lubrication_integral(x,z,n,t,h_coefficient_matrix,h0,l0)

% First we create the matrix SH_coeff
SH = zeros(n,n);
for l = 1: t-1
    SH(l,:) = x(l)^(3/2).*h_coefficient_matrix(3*l-2,:) ...
       + x(l)^(1/2).*h_coefficient_matrix(3*l-1,:) ...
       + h_coefficient_matrix(3*l,:);
end

for l = t:n
    SH(l,:) = x(l)^2.*h_coefficient_matrix(3*l-2,:) ...
       + x(l).*h_coefficient_matrix(3*l-1,:) ...
       + h_coefficient_matrix(3*l,:);
end

FSH = zeros(n,n);
for l = 1:n
    FSH(l,:) = l0*h0(l)^(-3) .* SH(l,:);
end

Q = zeros(n-1,n);

for k = 1:n-1
    Q(k,k) = (x(k+1)-z(k))^2 /(x(k+1)-x(k));
    Q(k,k+1) = (x(k+1)-z(k))*(z(k)-2*x(k)+x(k+1))/(x(k+1)-x(k));
    for r = k+1: n-1
        Q(k,r) = Q(k,r) - x(r);
        Q(k,r+1) = Q(k,r+1) + x(r+1);
    end
    Q(k,n) = Q(k,n) + (2/3)*x(n);
end

R = zeros(2*(n-1),2*n);
R(1:n-1,n+1:2*n) = Q*FSH;

end