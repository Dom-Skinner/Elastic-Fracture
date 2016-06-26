function deriv_pressure_num = pressure_map_derivative_num(x, z, h_coeffs,h_coefficient_matrix)

%returns the (n-1)xn matrix for the derivative of the map that sounds h' to
%pressure

%we first vary the 3n coefficients of h, to get a (n-1)x3n matrix. This is
%then multiplied by the map that takes h' to coefficients of h.
%h_coeffs = h_coefficient_matrix*hprime_old(n+1:2*n);
n = length(x);
t = round(n/2);
infinity = 10^(10);

%the coefficients of h
a = h_coeffs(1:3:3*n-2);
b = h_coeffs(2:3:3*n-1);
c = h_coeffs(3:3:3*n);

%functions to find (-2/h^3).


%we first vary a. varying the first t-1 values of a variation require
%integrating -2/h^3 against x^(3/2) in panel t. The next require
%integrating -2/h^3 against x^2.

dh_a_sqrt_fn = @(s,a,b,c) -2/(a*x^(3/2)+b*x^(1/2)+c)^3*x^(3/2);

dh_a_quad_fn = @(s,a,b,c) -2/(a*s^2+b*s+c)^3*s^2;

dh_b_sqrt_fn = @(s,a,b,c) -2/(a*s^(3/2)+b*s^(1/2)+c)^3*s^(1/2);

dh_b_quad_fn = @(s,a,b,c) -2/(a*s^2+b*s+c)^3*s;

dh_c_sqrt_fn = @(s,a,b,c) -2/(a*s^(3/2)+b*s^(1/2)+c)^3;

dh_c_quad_fn = @(s,a,b,c) -2/(a*s^2+b*s+c)^3;

%now produces a (n-1)x3n matrix with all these as entries

dp_num_h = zeros(n-1,3*n);

for j = 1:n-1
   for i = j:n
       %does the variation integral when i=j, integrate between z(j) and
       %x(i+1)
       if i == j
           %covers the case when i <= t-1
           if i <= t-1
               a_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3.*s.^(3/2);
               b_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3.*s.^(1/2);
               c_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3;
               dp_num_h(j,3*(i-1)+1) = integral(a_fn, z(j), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+2) = integral(b_fn, z(j), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+3) = integral(c_fn, z(j), x(i+1),'AbsTol',10^(-10));
           else
               a_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3.*s.^2;
               b_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3.*s.^1;
               c_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3;               
               dp_num_h(j,3*(i-1)+1) = integral(a_fn, z(j), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+2) = integral(b_fn, z(j), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+3) = integral(c_fn, z(j), x(i+1),'AbsTol',10^(-10));           
           end
       else
           if i <= t-1
               a_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3.*s.^(3/2);
               b_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3.*s.^(1/2);
               c_fn = @(s) -2./(a(i).*s.^(3/2)+b(i).*s.^(1/2)+c(i)).^3;               
               dp_num_h(j,3*(i-1)+1) = integral(a_fn, x(i), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+2) = integral(b_fn, x(i), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+3) = integral(c_fn, x(i), x(i+1),'AbsTol',10^(-10));
           elseif i ~= n
               a_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3.*s.^2;
               b_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3.*s.^1;
               c_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3;               
               dp_num_h(j,3*(i-1)+1) = integral(a_fn, x(i), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+2) = integral(b_fn, x(i), x(i+1),'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+3) = integral(c_fn, x(i), x(i+1),'AbsTol',10^(-10));            
           else
               a_fn = @(s) -2./((a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3).*s.^2;
               b_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3.*s.^1;
               c_fn = @(s) -2./(a(i).*s.^(2)+b(i).*s.^(1)+c(i)).^3;              
               dp_num_h(j,3*(i-1)+1) = integral(a_fn, x(i), infinity,'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+2) = integral(b_fn, x(i), infinity,'AbsTol',10^(-10));
               dp_num_h(j,3*(i-1)+3) = integral(c_fn, x(i), infinity,'AbsTol',10^(-10));                           
           end           
       end
   end
end

dp_num_h = real(dp_num_h);

deriv_pressure_num = dp_num_h; %*h_coefficient_matrix;