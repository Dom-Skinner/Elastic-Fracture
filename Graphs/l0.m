% This script is to find the parameters for the relationship
% $\lambda \approx \lambda_0 + D K_I^u $ where $\lambda_0,D$ are found 
% here numerically. 

%n_val = [350,407,465,524,815];
%x_val = [873,822,819,846,846];
clear
nstr = '524';
xstr = '846';
file = strcat('n',nstr,'x',xstr,'-mod');
load(file)
%load n400x40


s = 0.138673;
u = 4 - 6*s;
K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
p2 = polyfit(K(end-2:end).^u,lambda(end-2:end),2);

er    = abs(p1(2)-p2(3));
pc_er = abs(p1(2)-p2(3))/p1(2);

fprintf('Approximate interpolation error = %.2e\n',er)
fprintf('Approximate interpolation percentage error = %.2e%%\n',100*pc_er)

figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
lx = [K.^u 0];
hold on
plot(lx,p1(2)+p1(1).*lx,'k--', lx,p2(3)+p2(2).*lx+p2(1).*lx.^2,'b:', ...
    K.^u,lambda,'ko','LineWidth',2.5, 'MarkerFaceColor', 'k', 'MarkerSize', 10);

axis( [0, 0.5, 0.0555,0.0595]);
axis square
xlabel('$K_I^u$','Interpreter','Latex','fontsize',25)
ylabel('$\lambda$','Rotation',0, 'Position', [-0.08, 0.0574],...
    'Interpreter','Latex','fontsize',25)
set(gca,'fontsize',20')
set(gca,'TickLabelInterpreter', 'latex');
title('Plot of $K_I^u$ against $\lambda$',...
    'fontsize', 25,'Interpreter','latex');

s1 = strcat('Linear fit: $\lambda_0 =$',num2str(p1(2),5));
s2 = strcat('Quadratic fit: $\lambda_0 =$',num2str(p2(3),5));
text(0.27,0.0585,s1,'Interpreter','latex','fontsize',20)
text(0.27,0.0583,s2,'Interpreter','latex','fontsize',20)
text(0.27,0.0581,['$n=', nstr, '$'],'Interpreter','latex','fontsize',20)
text(0.27,0.0579,['$x_{end}=', xstr, '$'],'Interpreter','latex','fontsize',20)


legend({'Linear fit','Quadratic fit','$\lambda$'},'Interpreter','Latex'...
    ,'fontsize',20)
%return
fid = fopen('l0.csv','w');
fprintf(fid,'type,   lambda,    K \n');
for j = 1:numel(K)
    fprintf(fid, 'Kval,    %.5e,   %.5e \n',lambda(j), K(j).^u);
end
for j = 1:numel(lx)
    fprintf(fid, 'linear,    %.5e,   %.5e \n',(p1(2)+p1(1).*lx(j)),lx(j));
end
for j = 1:numel(lx)
    fprintf(fid, 'quadratic,    %.5e,   %.5e \n',(p2(3)+p2(2).*lx(j)+p2(1).*lx(j).^2),lx(j));
end

fclose(fid);

clear