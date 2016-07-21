clear plot
figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
ax = gca;
hold on
axis square

xlabel(ax,'$ \lambda $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ K_I $','Interpreter','latex','fontsize',25);
title(ax,'Scatter plot of toughness against speed with $n =800$, $x_{end}=30$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.006, 1])

clear data
load n800x30.mat
K = 3*sqrt(2*pi).*KI;
plot(ax,lambda,K, 'ko-')

export_fig ('K-lambda', '-pdf', '-transparent')
clear