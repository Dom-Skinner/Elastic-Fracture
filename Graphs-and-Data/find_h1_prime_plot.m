clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;

[h0_prime,~] = interpolate_hprime(x,n,hprime_data,K);
h1_prime = find_h1_prime(n,hprime_data,K,h0_prime);
h1_prime_data = zeros(n,numel(K));

for k = 1:n
    h1_prime_data(k,:) = (hprime_data(n+k,:)-h0_prime(k)).*K.^(-u);
end

figure('units','normalized','outerposition',[0 0 1 1]) % Makes figure fill 
% the whole screen. 
%


ax = gca;

dat = 1:199;
xplot(x(dat),x(dat).^(0.5-s).*h1_prime_data(dat,5)','o-',...
     x(dat),x(dat).^(0.5-s).*h1_prime_data(dat,6)','o-',...
     x(dat),x(dat).^(0.5-s).*h1_prime_data(dat,7)','o-',...
     x(dat),x(dat).^(0.5-s).*h1_prime_data(dat,9)','o-',...
     x(dat),x(dat).^(0.5-s).*h1_prime(dat),':','LineWidth',1.5,'MarkerSize', 3)
    
%axis([0, 0.15, 0.002, 0.0045])
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ h''x^{1/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $(h''-h_0'')x^{1-s}K_I^{-u}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.0014, 0.365])

s = '$K_I$ =';
s1 = strcat(s, num2str(K(5)));
s2 = strcat(s, num2str(K(6)));
s3 = strcat(s, num2str(K(7)));
s4 = strcat(s, num2str(K(9)));
legend({s1,s2,s3,s4,'Interpolated'},'Interpreter','latex')

clear