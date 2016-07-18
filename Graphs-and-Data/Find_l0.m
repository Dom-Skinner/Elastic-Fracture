% This script is to find the parameters for the relationship
% $\lambda \approx \lambda_0 + D K_I^u $ where $\lambda_0,D$ are found 
% here numerically. 

clear
load n800x30-extended.mat

s = 0.138673;
u = 4 - 6*s;
Kap = 3*sqrt(2*pi)*KI;
p = polyfit(Kap(end-13:end-3).^u,lambda(end-13:end-3),1);
D = p(1);
l0 = p(2);

figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%

plot(Kap.^u,l0+D.*Kap.^u, Kap.^u,lambda,'o--');
axis( [0, 1.5, 0.05,0.06]);
axis square
xlabel('$K_I^u$','Interpreter','Latex','fontsize',25)
set(gca,'fontsize',20')
set(gca,'TickLabelInterpreter', 'latex');

title('Plot of $K_I^u$ against $\lambda$',...
    'fontsize', 25,'Interpreter','latex');

legend({'$\lambda_0+DK_I^u$', '$\lambda$'},'Interpreter','Latex','fontsize',20)
fprintf('lambda_0 = %.2e, D = %.2e\n',l0,D);
fprintf('Where lambda = lambda_0 + D*K_I^u\n');

for l = 1:numel(lambda)-1
    p = polyfit(Kap(l:end).^u,lambda(l:end),1);
    D = p(1);
    l0 = p(2);
    fprintf('lambda_0 = %.4e, D = %.2e\n',l0,D);
end
