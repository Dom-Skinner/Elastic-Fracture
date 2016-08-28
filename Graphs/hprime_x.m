figure('units','normalized','outerposition',[0 0 0.5 1])
ax = gca;

load K_lambda
dat = 1:75;
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0 = p1(2);

[~,interp2] = interpolate_hprime(x,n,hprime_data,K,1/2,l0);
% Need to add /full fluid problem to path for this bit.

plot(ax,x(dat),hprime_data(n+dat,13)'.*x(dat).^(-1/6),'*-', ...
        x(dat),hprime_data(n+dat,14)'.*x(dat).^(-1/6),'*-', ...
        x(dat),hprime_data(n+dat,15)'.*x(dat).^(-1/6),'*-', ...
        x(dat),hprime_data(n+dat,17)'.*x(dat).^(-1/6),'*-', ...
        x(dat),interp2(dat).*x(dat).^(-1/6),'--',...
        'LineWidth',1.5,'MarkerSize', 8)

axis( [0.0, 0.015, 0.351,0.37]);
axis square
xlabel(ax,'$ x $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ h''x^{1/3} $','Interpreter','latex','fontsize',25);
title(ax,'Plot of $x$ against $h''x^{1/3}$',...
    'fontsize', 25,'Interpreter','latex');
set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')


s = '$K_I$ =';
s1 = strcat(s, num2str(K(13),3));
s2 = strcat(s, num2str(K(14),3));
s3 = strcat(s, num2str(K(15),3));
s4 = strcat(s, num2str(K(17),3));



legend({s1,s2,s3,s4,'Interpolated'}...
    ,'Interpreter','latex')


num = [{'one'}, {'two'}, {'three'}, {'four'}, {'five'}];
%%{
fid = fopen('hprime-x.csv','w');
fprintf(fid,'K,  x,   hprime,  \n');
for l = [1,2,3,5]
for j = dat
    fprintf(fid, [num{l},',    %.5e,    %.5e \n'],...
        x(j), hprime_data(n+j,12+l)*x(j).^(-1/6));
end
end

for j = dat
    fprintf(fid, 'interp,    %.5e,    %.5e \n',...
        x(j), interp2(j)*x(j).^(-1/6));
end

fclose(fid);
%}