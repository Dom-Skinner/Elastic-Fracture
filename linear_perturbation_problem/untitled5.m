clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.1386;
u = 4 - 6*s;
l0 = 0.05943;
D = -0.00809;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;


[h0_prime_l,~] = interpolate_hprime(x,n,hprime_data,K);

h_coefficient_matrix_l = hprime_to_h_l(x);
h_coefficient_matrix_s = hprime_to_h_s(x,s);

h0 = h_integrate(h0_prime_l',x,n,t,h_coefficient_matrix_l);



S_l = zeros(n,3*n);
S_s = zeros(n,3*n);
for k = 1:t-1
    S_l(3*k-2)=x(k)^(3/2);
    S_s(3*k-2)=x(k)^(s+1);
    
    S_l(3*k-1)=x(k)^(1/2);
    S_s(3*k-1)=x(k)^(s);
    
    S_l(3*k)=1;
    S_s(3*k)=1;
end

for k = t:n
    S_l(3*k-2)=x(k)^(2);
    S_s(3*k-2)=x(k)^(2);
    
    S_l(3*k-1)=x(k);
    S_s(3*k-1)=x(k);
    
    S_l(3*k)=1;
    S_s(3*k)=1;
end

SH_l = S_l*h_coefficient_matrix_l;
SH_s = S_s*h_coefficient_matrix_s;

vec =ones(1,n);
vec(1:t-1) = x(1:t-1);

I_l = diag(vec.^(0.5));
I_s = diag(vec.^(s));

Er = SH_l*I_l- SH_s*I_s;