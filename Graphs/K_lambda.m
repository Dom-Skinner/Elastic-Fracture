clear
load K_lambda

figure('units','normalized','outerposition',[0 0 0.5 1]) % Makes figure fill 
% the whole screen. Needed when using export fig.
%
ax = gca;
hold on
axis square

xlabel(ax,'$ \lambda $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ K_I $','Interpreter','latex','fontsize',25);
title(ax,['Scatter plot of toughness against speed with $n =', ...
    num2str(n),'$, $x_{end}=',num2str(round(x(end))),'$'],...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.006, 1])

plot(ax,lambda,K, 'ko-')

return
fid = fopen('K-lambda.csv','w');
fprintf(fid,'K,    lambda \n');
for j = 1:numel(K)
    fprintf(fid, '%.5e,   %.5e \n',K(j), lambda(j));
end
fclose(fid);
%export_fig ('K-lambda', '-pdf', '-transparent')
%clear