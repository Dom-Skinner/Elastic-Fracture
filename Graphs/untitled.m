clear
n_val = [200,400];
xend_val = [20,30,40,50,60,70];
for q  = 1:numel(n_val)
for l = 1:numel(xend_val)

file = strcat('n',num2str(n_val(q)),'x',num2str(xend_val(l)));
load(file)

u = 4 - 6*s;
K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0 = p1(2);
D_save(q,l) = p1(1);


er = ones(1,round(0.2*n));
for k = 20:round(0.25*n)
    p3 = polyfit(x(k:k+3) , H_LEFM_23(k:k+3)'.*x(k:k+3).^(-s),1);
    p4 = polyfit(x(k:k+4) ,H_LEFM_23(k:k+4)'.*x(k:k+4).^(-s) ,2);
    er(k) = abs(p3(2)-p4(3));
end
[~,I] = min(er(:));

p3 = polyfit(x(I:I+1) ,H_LEFM_23(I:I+1)'.*x(I:I+1).^(-s), 1);
l0_save(q,l) = l0;
l1_save(q,l) = -3*l0/p3(2);
clearvars -except l0_save l1_save D_save xend_val nstr l q n_val
end
end

figure('units','normalized','outerposition',[0 0 1 1])


ax1 = subplot(2,2,1);
plot(xend_val,l0_save(1,:),'o-',xend_val,l0_save(2,:),'o-')
xlabel(ax1,'$ x_{end} $','Interpreter','latex','fontsize',25);
ylabel(ax1,'$ \lambda_0$','Rotation',0,'Position', [10, 0.0594], ...
    'Interpreter','latex','fontsize',25);
title(ax1,'Numerical estimates of $\lambda_0$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax1,'TickLabelInterpreter', 'latex');
set(ax1,'fontsize',20')

ax2 = subplot(2,2,2);
plot(xend_val,l1_save(1,:),'o-',xend_val,l1_save(2,:),'o-')
xlabel(ax2,'$ x_{end} $','Interpreter','latex','fontsize',25);
ylabel(ax2,'$ \lambda_1$','Rotation',0,'Interpreter','latex','fontsize',25);
title(ax2,'Numerical estimates of $\lambda_1$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax2,'TickLabelInterpreter', 'latex');
set(ax2,'fontsize',20')

ax3 = subplot(2,2,3);
plot(xend_val,D_save(1,:),'o-',xend_val,D_save(2,:),'o-')
xlabel(ax3,'$ x_{end} $','Interpreter','latex','fontsize',25);
ylabel(ax3,'$ D$','Rotation',0,'Interpreter','latex','fontsize',25);
title(ax3,'Numerical estimates of $D$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax3,'TickLabelInterpreter', 'latex');
set(ax3,'fontsize',20')

ax4 = subplot(2,2,4);
s = 0.138673;
C = D_save./l1_save .*(l0_save).^(1-2*s);
plot(xend_val,C(1,:),'o-',xend_val,C(2,:),'o-')
xlabel(ax4,'$ x_{end} $','Interpreter','latex','fontsize',25);
ylabel(ax4,'$C$','Rotation',0,'Interpreter','latex','fontsize',25);
title(ax4,'Numerical estimates of $C$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax4,'TickLabelInterpreter', 'latex');
set(ax4,'fontsize',20')
%fprintf('\n\n xend = %.3e%%\n',xend)
%fprintf('l0 = %.3e%%\n',l0)
%fprintf('l1 = %.3e%%\n',-3*l0/p3(2))
