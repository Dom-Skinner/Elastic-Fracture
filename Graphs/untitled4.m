% This script is a check on the lubrication matrix R.
% We can find \tilde{\Pi}' via interpolating from p. Alternatively, we can
% find \tilde{P} by first finding \tilde{H}, and then integrating the 
% lubrication equation, which is equivalent to multiplying by the matrix R
% If both give the same answer, that is good news.


clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;
l0 = 0.05943;
D = -0.00809;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;


[~,h0_prime] = interpolate_hprime(x,n,hprime_data,K);
h_coefficient_matrix = hprime_to_h_s(x,0.5);
h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5);
h0_z = h_integrate(h0_prime',z,n-1,t,h_coefficient_matrix,0.5);
h1_prime = find_h1_prime(n,hprime_data,K,h0_prime);


pprime_data=find_pprime(lambda,x,hprime_data,n,t);
p0_prime = find_p0_prime(n,pprime_data,K);
p1_prime = find_p1_prime(n,pprime_data,K,p0_prime);


R2 = lubrication_integral(x,z,n,t,h_coefficient_matrix,h0,h0_z,l0,0.5);

R = lubrication_integral_s(x,z,n,t,h_coefficient_matrix,h0,l0,0.5);
R=R(1:n-1,n+1:2*n);


Hprime_tilde = h0_prime - (3*l0/D)*h1_prime;
Pprime_tilde = p0_prime - (3*l0/D)*p1_prime;

P_s = R2*Hprime_tilde';
Pprime_s = (P_s(2:end)-P_s(1:end-1))'./(z(2:end)-z(1:end-1));

P = R*Hprime_tilde';
Pprime = (P(2:end)-P(1:end-1))'./(z(2:end)-z(1:end-1));

plot(z(1:n-2),z(1:n-2).^(2-s).*Pprime_s,'o',z(1:n-2),z(1:n-2).^(2-s).*Pprime,'o',...
    x,x.^(2-s).*Pprime_tilde)

