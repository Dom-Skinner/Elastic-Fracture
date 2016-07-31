clear plot
figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. 
%
ax = gca;
hold on

clear data
load n200x30-tmp %n400x50-extended.mat
dat = 1:80;
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);

[interp1,interp2] = interpolate_hprime(x,n,hprime_data,K,1/2,p1(2));

plot(ax,x(dat),hprime_data(n+dat,5)'.*x(dat).^(-1/6), '*-', ...
    x(dat),hprime_data(n+dat,6)'.*x(dat).^(-1/6),'*-', ...
            x(dat),hprime_data(n+dat,7)'.*x(dat).^(-1/6),'*-', ...
x(dat),interp1(dat).*x(dat).^(-1/6),':',...
x(dat),interp2(dat).*x(dat).^(-1/6),'--',...
'LineWidth',1.5,'MarkerSize', 8)

axis( [0.0, 0.01, 0.35,0.38]);
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ h''x^{1/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $h''x^{1/3}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.0014, 0.365])

s = '$K_I$ =';
s1 = strcat(s, num2str(K(5)));
s2 = strcat(s, num2str(K(6)));
s3 = strcat(s, num2str(K(7)));


legend({s1,s2,s3,'Interpolated','Interpolated with LEFM correction'}...
    ,'Interpreter','latex')

clear