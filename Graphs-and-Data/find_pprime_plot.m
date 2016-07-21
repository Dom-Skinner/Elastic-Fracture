clear plot
figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. 
%
ax = gca;
hold on

clear data
load n400x50-extended.mat
dat = 1:80;
K_val=[1, 4, 8];
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
pprime_data=find_pprime(lambda,x,hprime_data,n,t);

p0_prime = find_p0_prime(n,pprime_data,K);

plot(ax,x(dat),pprime_data(dat,K_val(1))'.*x(dat).^(4/3), '*-', ...
        x(dat),pprime_data(dat,K_val(2))'.*x(dat).^(4/3),'*-', ...
        x(dat),pprime_data(dat,K_val(3))'.*x(dat).^(4/3) ,'*-', ...
        x(dat),p0_prime(dat).*x(dat).^(4/3),':',...
        'LineWidth',1.5,'MarkerSize', 8)

axis( [0.0, 0.06, 0,0.25]);
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ p''x^{4/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $p''x^{4/3}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.006, 0.16])

s = '$K_I$ =';
s1 = strcat(s, num2str(K(K_val(1))));
s2 = strcat(s, num2str(K(K_val(2))));
s3 = strcat(s, num2str(K(K_val(3))));


legend({s1,s2,s3,'Interpolated'}...
    ,'Interpreter','latex')

clear