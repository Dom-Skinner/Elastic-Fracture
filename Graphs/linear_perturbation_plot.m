clear
load n1022x773
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0 = 0.0591;
%{
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

[h0_prime ,h0_prime_LEFM] = interpolate_hprime(x,n,hprime_data,K,0.5,l0);

h_coefficient_matrix = hprime_to_h_s(x,0.5);

h0 = h_integrate(h0_prime',x,n,t,h_coefficient_matrix,0.5 );
h0_z = h_integrate(h0_prime',z,n-1,t,h_coefficient_matrix,0.5 );

h0_LEFM = h_integrate(h0_prime_LEFM',x,n,t,h_coefficient_matrix,0.5 );
h0_LEFM_z = h_integrate(h0_prime_LEFM',z,n-1,t,h_coefficient_matrix,0.5 );



h0_prime_LEFM_23 = convert(0.5,2/3,n,t,x,h0_prime_LEFM);
h_coefficient_matrix_23 = hprime_to_h_s(x,2/3);
h0_LEFM_23 = h_integrate(h0_prime_LEFM_23,x,n,t, ... 
    h_coefficient_matrix_23,2/3);
h0_LEFM_23_z = h_integrate(h0_prime_LEFM_23,z,n-1,t, ...
    h_coefficient_matrix_23,2/3);

[~,H] = linear_perturbation_solve(n,t,xmax, h0,h0_z,l0,kernel_matrix_s,s);

[~,H_LEFM] = linear_perturbation_solve(n,t,xmax, h0_LEFM,h0_LEFM_z,l0,...
    kernel_matrix_s,s);
%}
[~,H_LEFM_23] = linear_perturbation_solve(n,t,x,z, h0_LEFM_23,...
    h0_LEFM_23_z,l0,kernel_matrix,s);


%figure('units','normalized','outerposition',[0 0 0.5 1])
%plot(x,H'.*x.^-s,'o',x,H_LEFM'.*x.^-s,'o',x,H_LEFM_23'.*x.^-s,'*')
hold on
plot(x,H_LEFM_23'.*x.^-s,'*')


ax = gca;
axis square
xlabel(ax,'$ \xi $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \tilde{H} \xi^{-s}$','Interpreter','latex','fontsize',25);
title(ax,'Linear perturbation problem estimates of $\tilde{H}\xi^{-s}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
legend({'$H_0$ interpolated','$H_0$ with LEFM correction',...
    '$H_0$ with LEFM correction and $2/3$ integration'...
    },'Interpreter','Latex'...
    ,'fontsize',20,'Location','southeast')

axis([0,0.15,0.5,0.75])
%export_fig ('linear-perturbation', '-pdf', '-transparent')

er = ones(1,round(0.2*n));
for k = round(0.05*n):round(0.2*n)
    p3 = polyfit(x(k:k+3) , H_LEFM_23(k:k+3)'.*x(k:k+3).^(-s),1);
    p4 = polyfit(x(k:k+4) ,H_LEFM_23(k:k+4)'.*x(k:k+4).^(-s) ,2);
    er(k) = abs(p3(2)-p4(3));
end
[~,I] = min(er(:));

p3 = polyfit(x(I:I+3) ,H_LEFM_23(I:I+3)'.*x(I:I+3).^(-s), 1);
p4 = polyfit(x(I:I+4) ,H_LEFM_23(I:I+4)'.*x(I:I+4).^(-s), 2);

%interp2(1:I) = (p3(1)*x(1:I)+p3(2)).*x(1:I).^(2/3-a);
fprintf('intercept = %.2e%%\n',p3(2))



