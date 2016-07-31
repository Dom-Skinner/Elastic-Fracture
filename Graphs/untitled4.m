clear
n_val = [300,400];
xend_val = [20,24,32,38,45,48];

for k1 = 1:2
for l = 1:numel(xend_val)
 n = n_val(k1);
file = strcat('n',num2str(n),'x',num2str(xend_val(l)),'-tmp');
load(file)

    s = 0.138673;

u = 4 - 6*s;
K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0 = p1(2);
D_save(l,k1) = p1(1);


l0_save(l,k1) = l0;
clearvars -except l0_save D_save xend_val n_val nstr l  n k1
end
end

figure('units','normalized','outerposition',[0 0 1 1])

hold on
plot(xend_val.^(-2),l0_save(:,1),'o-',xend_val.^(-2),l0_save(:,2),'o-')

xlabel('$ 1/x_{end}^2 $','Interpreter','latex','fontsize',25);
ylabel('$ \lambda_0$','Rotation',0,'Position', [10, 0.0594], ...
    'Interpreter','latex','fontsize',25);
title('Numerical estimates of $\lambda_0$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');

axis([0,0.0025,0.059,0.0595])