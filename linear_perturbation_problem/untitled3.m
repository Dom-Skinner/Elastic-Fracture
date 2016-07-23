clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;
l0 = 0.05943;
D = -0.00809;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

[h0_prime,~] = interpolate_hprime(x,n,hprime_data,K);


h_coefficient_matrix = hprime_to_h_s(x,0.5);
h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5);

h0_prime_s = convert(0.5,2/3,n,t,x,h0_prime);
h_coefficient_matrix_s = hprime_to_h_s(x,2/3);

h0_s = h_integrate(h0_prime_s,x,n,t,h_coefficient_matrix_s,2/3);

plot(x,h0,x,h0_s)
