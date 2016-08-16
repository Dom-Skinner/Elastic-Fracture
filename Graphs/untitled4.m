clear
n_val = [200,400,500,600,800];
%xend_val = [20,24,32,38,45,48];
xend_val = [20,30,40,50,60,70];


for k1 = 1:numel(n_val)
    for l = 1:numel(xend_val)
        n = n_val(k1);
        if n_val(k1) ~= 500 && n_val(k1) ~= 600
            file = strcat('n',num2str(n),'x',num2str(xend_val(l)));
        else
            file = strcat('n',num2str(n),'x',num2str(xend_val(l)),'-tmp');
        end
        load(file)

        s = 0.138673;

        u = 4 - 6*s;
        K = 3*sqrt(2*pi)*KI;
        p1 = polyfit(K(1:end).^u,lambda(1:end),1);
        l0 = p1(2);
        
        D_save(l,k1) = p1(1);
        l0_save(l,k1) = l0;
        
        clearvars -except l0_save D_save xend_val n_val nstr l  n k1
    end
end




figure('units','normalized','outerposition',[0 0 1 1])

hold on
%{
plot(xend_val.^(-1),D_save(:,1),'o-',xend_val.^(-1),D_save(:,2),'+-',...
    xend_val.^(-1),D_save(:,3),'*-',xend_val.^(-1),D_save(:,4),'^-',...
    xend_val.^(-1),D_save(:,5),'o-')
%}
%%{
plot(xend_val.^(-2),l0_save(:,1),'o-',xend_val.^(-2),l0_save(:,2),'o-',...
    xend_val.^(-2),l0_save(:,3),'o-',xend_val.^(-2),l0_save(:,4),'o-',...
    xend_val.^(-2),l0_save(:,5),'o-')
%}


xlabel('$ 1/x_{end} $','Interpreter','latex','fontsize',25);
ylabel('$ \lambda_0$','Rotation',0,'Position', [0.006, -0.0079], ...
    'Interpreter','latex','fontsize',25);
title('Numerical estimates of $\lambda_0$ against $x_{end}$',...
    'fontsize', 25,'Interpreter','latex');
legend({'n=200','n=400','n=500, (New)', 'n=600, (New)', 'n=800'},...
    'Interpreter','latex','fontsize',25)
%axis([0,0.0025,0.059,0.0597])