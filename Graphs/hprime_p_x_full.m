figure('units','normalized','outerposition',[0 0 0.5 1])

clear

load n500x50-modified
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
for l = 1:numel(lambda)
h_coefficient_matrix = hprime_to_h(x,t,lambda(l));
h_coeffs = h_coefficient_matrix*hprime_data(n+1:2*n,l);
%finds the pressure: an n-1 vector
pressure(:,l) = pprime_to_p(x,z,h_coeffs,t);
hprime_data(n+1:n+t-1,l) = x(1:t-1)'.^(-0.5).*hprime_data(n+1:n+t-1,l);
end



% Different K values are indistiguishable at this scale, can check this by
% running the following:
% plot(ax,x,hprime_data(n+1:end,1),x,hprime_data(n+1:end,4), ...
%            x,hprime_data(n+1:end,7))
subplot(1,2,1)
ax = gca;
plot(ax,x,hprime_data(n+1:end,1),'LineWidth',2)

axis( [0, 4, 0,5]);
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ h'' $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.5, 2.5])

subplot(1,2,2)
ax = gca;
plot(z,pressure(:,1),'LineWidth',2)
axis( [0, 4, -3,0]);
axis square

xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ p $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.4, -1.5])
