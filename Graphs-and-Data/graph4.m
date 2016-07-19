clear plot
figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
ax = gca;
hold on

clear data
load n400x40.mat
dat = 1:40;
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;

[interp1, interp2] = interpolate_hprime(x,n,hprime_data,K);

plot(ax,x(dat),hprime_data(n+dat,9)'.*x(dat).^(-1/6), '*-', ...
    x(dat),hprime_data(n+dat,10)'.*x(dat).^(-1/6),'*-', ...
            x(dat),hprime_data(n+dat,11)'.*x(dat).^(-1/6),'*-', ...
x(dat),interp1(dat).*x(dat).^(-1/6),x(dat),interp2(dat).*x(dat).^(-1/6))

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
s1 = strcat(s, num2str(K(9)));
s2 = strcat(s, num2str(K(10)));
s3 = strcat(s, num2str(K(11)));


legend({s1,s2,s3,'Interpolated','Interpolated with LEFM correction'}...
    ,'Interpreter','latex')

