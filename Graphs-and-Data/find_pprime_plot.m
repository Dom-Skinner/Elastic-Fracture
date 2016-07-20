clear plot
figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. 
%
ax = gca;
hold on

clear data
load n400x50-extended.mat
dat = 1:80;
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
pprime_data=find_pprime(lambda,x,hprime_data,n,t);

[interp1, interp2] = interpolate_hprime(x,n,hprime_data,K);

plot(ax,x(dat),pprime_data(dat,1)'.*x(dat).^(5/6), '*-', ...
    x(dat),pprime_data(dat,4)'.*x(dat).^(5/6),'*-', ...
            x(dat),pprime_data(dat,8)'.*x(dat).^(5/6) ,'*-')%, ...
%x(dat),interp1(dat).*x(dat).^(-1/6),':',...
%x(dat),interp2(dat).*x(dat).^(-1/6),'--',...
%'LineWidth',1.5,'MarkerSize', 8)

axis( [0.0, 0.06, 0,3]);
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ p''x^{4/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $p''x^{4/3}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.006, 1.5])

s = '$K_I$ =';
s1 = strcat(s, num2str(K(6)));
s2 = strcat(s, num2str(K(7)));
s3 = strcat(s, num2str(K(8)));


legend({s1,s2,s3,'Interpolated'}...
    ,'Interpreter','latex')

clear