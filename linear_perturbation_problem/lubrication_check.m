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


[h0_prime,~] = interpolate_hprime(x,n,hprime_data,K);
h_coefficient_matrix = hprime_to_h_s(x,0.5);
h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5);
h1_prime = find_h1_prime(n,hprime_data,K,h0_prime);
h1 = h_integrate(h1_prime',x,n,t,h_coefficient_matrix,0.5);


pprime_data=find_pprime(lambda,x,hprime_data,n,t);
p0_prime = find_p0_prime(n,pprime_data,K);
p1_prime = find_p1_prime(n,pprime_data,K,p0_prime);

h_coefficient_matrix_s = hprime_to_h_s(x,s);


R_s = lubrication_integral_s(x,z,n,t,h_coefficient_matrix_s,h0,l0,s);
R_s=R_s(1:n-1,n+1:2*n);

R = lubrication_integral_s(x,z,n,t,h_coefficient_matrix,h0,l0,0.5);
R=R(1:n-1,n+1:2*n);


Hprime_tilde = h0_prime - (3*l0/D)*h1_prime;
H_tilde = h0 - (3*l0/D)*h1;
Pprime_tilde = p0_prime - (3*l0/D)*p1_prime;

Hprime_tilde_s = convert(0.5,s,n,t,x,Hprime_tilde);
P_s = R_s*Hprime_tilde_s;
Pprime_s = (P_s(2:end)-P_s(1:end-1))'./(z(2:end)-z(1:end-1));

P = R*Hprime_tilde';
Pprime = (P(2:end)-P(1:end-1))'./(z(2:end)-z(1:end-1));

figure('units','normalized','outerposition',[0 0 0.5 1])
plot(z(1:n-2),z(1:n-2).^(2-s).*Pprime_s,'o',z(1:n-2),z(1:n-2).^(2-s).*Pprime,'o',...
    x,x.^(2-s).*Pprime_tilde)


ax = gca;
axis square
xlabel(ax,'$ \xi $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \tilde{\Pi}''\xi^{2-s} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of alternative ways to calculate $\tilde{\Pi}''$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
legend({'$x^s$ linearisation','$x^{1/2}$ linearisation'...
    'Via interpolation of $\Pi(\xi;\kappa)$'},'Interpreter','Latex'...
    ,'fontsize',20,'Location','southeast')


axis([0,30, -0.45,0.05])


