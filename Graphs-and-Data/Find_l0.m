% This script is to find the parameters for the relationship
% $\lambda \approx \lambda_0 + D K_I^u $ where $\lambda_0,D$ are found 
% here numerically. 

clear
load n400x40.mat

s = 0.138673;
u = 4 - 6*s;
Kap = 3*sqrt(2*pi)*KI;
p1 = polyfit(Kap.^u,lambda,1);
p2 = polyfit(Kap(end-7:end-3).^u,lambda(end-7:end-3),1);
p3 = polyfit(Kap(end-3:end).^u,lambda(end-3:end),1);


figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%

plot(Kap.^u,lambda,'o', Kap.^u,p1(2)+p1(1).*Kap.^u,  ...
Kap.^u,p2(2)+p2(1).*Kap.^u, Kap.^u,p3(2)+p3(1).*Kap.^u);

axis( [0, 1.5, 0.05,0.06]);
axis square
xlabel('$K_I^u$','Interpreter','Latex','fontsize',25)
set(gca,'fontsize',20')
set(gca,'TickLabelInterpreter', 'latex');
title('Plot of $K_I^u$ against $\lambda$',...
    'fontsize', 25,'Interpreter','latex');

s1 = strcat(num2str(p1(2),4),num2str(p1(1),4),'$K_I^u$');
s2 = strcat(num2str(p2(2),4),num2str(p2(1),4),'$K_I^u$');
s3 = strcat(num2str(p3(2),4),num2str(p3(1),4),'$K_I^u$');



legend({'$\lambda$',s1,s2,s3},'Interpreter','Latex','fontsize',20)


clear s1 s2 s3 p1 p2 p3 u s Kap