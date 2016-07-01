clear plot
figure('units','normalized','outerposition',[0 0 1 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%


for n=1:2
	ax(n) = subplot(1,2,n);
    hold(ax(n),'on');

end
%
clear data
load n800x30.mat
K = 3*sqrt(2*pi).*KI;
plot(ax(1),lambda,K.^(3.168),'black')
scatter(ax(1),lambda,K.^(3.168),70,'black')

plot(ax(2),lambda(8:14),K(8:14).^(3.168),'black')
scatter(ax(2),lambda(8:14),K(8:14).^(3.168),70,'black')

for n=1:2
	xlabel(ax(n),'$ \lambda $','Interpreter','latex','fontsize',25);
    ylabel(ax(n),'$ K_I^{3.16\dots} $','Interpreter','latex','fontsize',25);
	set(ax(n),'fontsize',20')
	pbaspect(ax(n), [1 1 1]); % So that the plots appear square
	set(ax(n),'TickLabelInterpreter', 'latex');
	set(ax(n), 'Color', 'none')

end
%
title(ax(1),'Plot of suspected relationship between $\lambda$ and $K_I$','fontsize',...
    25,'Interpreter','latex');
title(ax(2),'Zoom of plot, for $n=800$, $x_{max}=30$','fontsize',25,...
    'Interpreter','latex');
%
set(get(ax(1),'YLabel'),'Rotation',0,'Position', [-0.007, 4.5])
set(get(ax(2),'YLabel'),'Rotation',0,'Position', [0.0566, 0.15])

%

%title(ax,'Scatter plot of toughness against speed with $n =800$, $x_{end}=30$',...
%    'fontsize', 25,'Interpreter','latex');
%set(ax,'TickLabelInterpreter', 'latex');
%set(ax,'fontsize',20')
%set(get(ax,'YLabel'),'Rotation',0)
%set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.006, 1])



export_fig ('sus-relation', '-png', '-transparent','-m1.5')