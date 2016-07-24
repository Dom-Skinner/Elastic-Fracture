clear
load n400x40-extended

clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
l0 = 0.05943;
D = -0.00809;


[~,h0_prime] = interpolate_hprime(x,n,hprime_data,K);
h_coefficient_matrix = hprime_to_h_s(x,0.5);
h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5);

h1_prime = find_h1_prime(n,hprime_data,K,h0_prime);

h_coefficient_matrix_s = hprime_to_h_s(x,s);


Hprime_tilde = h0_prime - (3*l0/D)*h1_prime;

Hprime_tilde_s = convert(0.5,s,n,t,x,Hprime_tilde);


p1 = polyfit(x(25:35) , Hprime_tilde_s(25:35)',1);
Hprime_tilde_s(1:20)=p1(1).*x(1:20)+p1(2);



plot(x,Hprime_tilde_s,'o',x,Hprime_tilde_u)
axis([0,0.1,0.05,0.2])