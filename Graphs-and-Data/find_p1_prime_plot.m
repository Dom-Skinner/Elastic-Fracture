clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;

pprime_data=find_pprime(lambda,x,hprime_data,n,t);
p0_prime = find_p0_prime(n,pprime_data,K);

p1_prime = find_p1_prime(n,pprime_data,K,p0_prime);
p1_prime_data = zeros(n,numel(K));

for k = 1:n
    p1_prime_data(k,:) = (pprime_data(k,:)-p0_prime(k)).*K.^(-u);
end

figure('units','normalized','outerposition',[0 0 1 1]) % Makes figure fill 
% the whole screen. 
%

ax = gca;

dat = 2:199;
plot(x(dat),x(dat).^(2-s).*p1_prime_data(dat,5)','o-',...
     x(dat),x(dat).^(2-s).*p1_prime_data(dat,6)','o-',...
     x(dat),x(dat).^(2-s).*p1_prime_data(dat,7)','o-',...
     x(dat),x(dat).^(2-s).*p1_prime(dat),':','LineWidth',1.5,'MarkerSize', 3)
    
axis([0, 0.3, -0.021, -0.015])
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $(p''-p_0'')x^{1-s}K_I^{-u}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')

s = '$K_I$ =';
s1 = strcat(s, num2str(K(5)));
s2 = strcat(s, num2str(K(6)));
s3 = strcat(s, num2str(K(7)));
legend({s1,s2,s3,'Interpolated'},'Interpreter','latex')

clear