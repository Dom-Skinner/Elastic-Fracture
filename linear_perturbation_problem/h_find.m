% This script finds h0 given hprime_data.
% Saves the fancy, interpolated + LEFM + 2/3 integration estimate.
clear
file = 'n700x875';
load(file)
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
l0 = 0.0591;


[h0_prime ,h0_prime_LEFM] = interpolate_hprime(x,n,hprime_data,K,0.5,l0);

h_coefficient_matrix = hprime_to_h_s(x,0.5,t);

h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5 );
h0_z = h_integrate(h0_prime',z,n-1,t,h_coefficient_matrix,0.5 );

h0_LEFM = h_integrate(h0_prime_LEFM',x,n,t,h_coefficient_matrix,0.5 );
h0_LEFM_z = h_integrate(h0_prime_LEFM',z,n-1,t,h_coefficient_matrix,0.5 );


h0_prime_LEFM_23 = convert(0.5,2/3,n,t,x,h0_prime_LEFM);
h_coefficient_matrix_23 = hprime_to_h_s(x,2/3,t);
h0_LEFM_23 = h_integrate(h0_prime_LEFM_23,x,n,t, ... 
    h_coefficient_matrix_23,2/3);
h0_LEFM_23_z = h_integrate(h0_prime_LEFM_23,z,n-1,t, ...
    h_coefficient_matrix_23,2/3);

plot(x,h0_LEFM,'*-',x,h0,'+-',x,h0_LEFM_23,'o-');

%save(file,'-append','h0_LEFM_23','h0_LEFM_23_z')