clear plot
figure('units','normalized','outerposition',[0 0 1 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
ax = gca;
hold on

axis( [0.0, 0.3, 0.3,0.4]);
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ h''x^{1/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $h''x^{1/3}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.02, 0.35])

clear data
load n800x40.mat
dat = 1:300;
plot(ax,x(dat),hprime_data(dat,2)'.*x(dat).^(-1/6),...
    x(dat),hprime_data(dat,6)'.*x(dat).^(-1/6), ...
x(dat),hprime_data(dat,14)'.*x(dat).^(-1/6))
s = '$K_I$ =';
s1 = strcat(s, num2str(2*sqrt(3*pi)*KI(2)));
s2 = strcat(s, num2str(2*sqrt(3*pi)*KI(6)));
s3 = strcat(s, num2str(2*sqrt(3*pi)*KI(14)));

legend({s1,s2,s3},'Interpreter','latex')
    


export_fig ('hprime-x', '-png', '-transparent','-m1.5')