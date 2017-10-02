figure('units','normalized','outerposition',[0 0 1 1])

clear

load K_lambda_rescaled
s = 0.138673;
u = 4 - 6*s;

K = 3*sqrt(2*pi)*KI;
for l = 1:numel(lambda)
h_coefficient_matrix = hprime_to_h(x,t,lambda(l));
h_coeffs = h_coefficient_matrix*hprime_data(n+1:2*n,l);
%finds the pressure: an n-1 vector
p(:,l) = pprime_to_p(x,z,h_coeffs,t);
h(:,l) = h_integrate(hprime_data(n+1:2*n,l),x,n,t, ... 
    h_coefficient_matrix,1/2);

hprime_data(n+1:n+t-1,l) = x(1:t-1)'.^(-0.5).*hprime_data(n+1:n+t-1,l);
hprime_data(1:t-1,l) = x(1:t-1)'.^(-0.5).*hprime_data(1:t-1,l);


end

hprime = hprime_data(n+1:end,:);
gprime = hprime_data(1:n,:);

dat = [1:2:10,  11:4:40,45:10:numel(x)];
x = x(dat);
hprime = hprime(dat,:);
gprime = gprime(dat,:);
dat(end)=[];
z = z(dat);
p = p(dat,:);


% Different K values are indistiguishable at this scale, can check this by
% running the following:
% plot(ax,x,hprime_data(n+1:end,1),x,hprime_data(n+1:end,4), ...
%            x,hprime_data(n+1:end,7))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
ax = gca;

hold on
p1 = plot(ax,x,hprime(:,2),'LineWidth',2);
p2 = plot(x,hprime(:,10),':','LineWidth',2);

axis( [0, 2, 0.8,2.5]);
axis square
xlabel(ax,'$ \xi $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ H'' $','Interpreter','latex','fontsize',25);
legend([p1,p2],['\kappa_I =', num2str(K(2),3)],...
    ['\kappa_I =',num2str(K(10),2)],'Location','SouthEast')

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
%set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.5, 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4)
ax = gca;
k1 = log(KI(1));
k2 = log(KI(10));

plot(ax,[0,6],[0,6],'k',[-20,-8],[10+k1,4+k1],'k',[-20,-8],[10+k2,4+k2],'k',...
'LineWidth',1.5)
hold on
p1 = plot(ax,log(x),log(hprime(:,2)),'+','MarkerFaceColor','b');
p2 = plot(log(x),log(hprime(:,10)),'d','LineWidth',2);

axis( [-15.5, 6, -0.4,6.5]);
axis square
xlabel(ax,'$ \log(\xi) $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \log(H'') $','Interpreter','latex','fontsize',25);
legend([p1,p2],['\kappa_I =', num2str(K(2),3)], ['\kappa_I =',num2str(K(10),2)])

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
%set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.5, 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2)
ax = gca;
plot(ax,x,gprime(:,1),x,gprime(:,10),':','LineWidth',2)
axis( [0, 2, 0.5,1.5]);
axis square

xlabel(ax,'$ \xi $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ G'' $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.4, 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,5)
ax = gca;
plot(ax,log(x),log(gprime(:,1)),'+',log(x),log(gprime(:,10)),'d','LineWidth',2)
hold on
c = 0.4;
plot(ax,[-8, -2], [2+c,-1+c],'k','LineWidth',1.5)
axis( [-6, 3, -0.8,1]);
axis square

xlabel(ax,'$ \log(\xi) $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \log(G'') $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
%set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.4, 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,3)
ax = gca;

plot(z,lambda(2)*p(:,2),z,lambda(10)*p(:,10),':','LineWidth',2)
axis( [0, 2, -1,0]);
axis square

xlabel(ax,'$ \xi $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \Pi $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.4, -1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6)
ax = gca;

plot(ax,log(z),lambda(2)*p(:,2),'+', ...
     log(z),lambda(10)*p(:,10),'d',...
     [-16,-6],lambda(2)*9*pi/(2*K(2)^2)*[-16,-6], 'k', ...
    [-16,-8],lambda(10)*9*pi/(2*K(10)^2)*[-16,-8]+ 7.4, 'k','LineWidth',2)
axis( [-16, 1, -16,1]);
axis square

xlabel(ax,'$ \log(\xi) $','Interpreter','latex','fontsize',25);
ylabel(ax,'$ \Pi $','Interpreter','latex','fontsize',25);

set(ax,'TickLabelInterpreter', 'latex');
set(ax,'fontsize',20')
set(get(ax,'YLabel'),'Rotation',0)
%set(get(ax,'YLabel'),'Rotation',0, 'Position', [-0.4, -1.5])


fid = fopen('hprime-p-x-full-1.csv','w');
fprintf(fid,'K,  x,   hprime,  gprime,   \n');
for l = [2,10]
for j = 1:numel(x)
    fprintf(fid, '%.5e,    %.5e,   %.5e,    %.5e,  \n',...
        K(l),x(j), hprime(j,l), gprime(j,l) );
end
end
fclose(fid);


fid = fopen('hprime-p-x-full-2.csv','w');
fprintf(fid,'K,  z,   pressure\n');
for l = [2,10]
for j = 1:numel(z)
    fprintf(fid, '%.5e,    %.5e,   %.5e \n',...
        K(l),z(j), lambda(l)*p(j,l) );
end
end
fclose(fid);
