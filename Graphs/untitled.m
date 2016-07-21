clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;

[h0_prime,~] = interpolate_hprime(x,n,hprime_data,K);
h_coefficient_matrix = hprime_to_h_l(x);
h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix);
h1_prime = find_h1_prime(n,hprime_data,K,h0_prime);
h1 = h_integrate(h1_prime',x,n,t,h_coefficient_matrix);

pprime_data=find_pprime(lambda,x,hprime_data,n,t);
p0_prime = find_p0_prime(n,pprime_data,K);

p1_prime = find_p1_prime(n,pprime_data,K,p0_prime);

% note that neither p1_prime nor p0_prime, nor h0,h1 are stored with
% the built in x^(-0.5) singularity

dat = 1:400;
plot(x(dat),h0(dat)'.^2.*p1_prime(dat)+2*h0(dat)'.*h1(dat)'.*p0_prime(dat))