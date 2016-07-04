clear plot
figure('units','normalized','outerposition',[0 0 1 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
ax = gca;
hold on

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

clear data
load n400x40.mat
dat = 1:30;
plot(ax,x(dat),hprime_data(401:430,9)'.*x(dat).^(-1/6), '*-', ...
    x(dat),hprime_data(400+dat,10)'.*x(dat).^(-1/6),'*-', ...
            x(dat),hprime_data(400+dat,11)'.*x(dat).^(-1/6),'*-', ...
x(dat),hprime_data(400+dat,12)'.*(x(dat).^(-1/6)),'*-')
s = '$K_I$ =';
s1 = strcat(s, num2str(3*sqrt(2*pi)*KI(9)));
s2 = strcat(s, num2str(3*sqrt(2*pi)*KI(10)));
s3 = strcat(s, num2str(3*sqrt(2*pi)*KI(11)));
s4 = strcat(s, num2str(3*sqrt(2*pi)*KI(12)));

legend({s1,s2,s3,s4},'Interpreter','latex')
    


%export_fig ('hprime-x', '-png', '-transparent','-m1.5')