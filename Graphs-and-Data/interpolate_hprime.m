%clear data
%load n400x40.mat

%K = 3*sqrt(2*pi)*KI;
function [interp1, interp2] = interpolate_hprime(x,n,hprime_data,K,fit)
s = 0.138673;
u = 4 - 6*s;
interp1 = zeros(1,n);
er = zeros(1,n);
pc_er = zeros(1,n);

% First go at interpolation of hprime
for l = 1:n
    p1 = polyfit(K(fit:fit+1).^u , hprime_data(n+l,fit:fit+1),1);
    p2 = polyfit(K(fit:fit+2).^u , hprime_data(n+l,fit:fit+2),2);
    interp1(l) = p1(2);
    er(l)    = abs(p1(2)-p2(3));
    pc_er(l) = abs(p1(2)-p2(3))/p1(2);
end

fprintf('Approximate interpolation error = %.2e%%\n',100*max(pc_er(40:end)))
% Now we correct for the LEFM boundary
p3 = polyfit(x(20:21) , interp1(20:21).*x(20:21).^(-1/6),1);
p4 = polyfit(x(20:22) , interp1(20:22).*x(20:22).^(-1/6),2);
interp2 = interp1;
interp2(1:20) = (p3(1)*x(1:20)+p3(2)).*x(1:20).^(1/6);
fprintf('Approximate LEFM correction error = %.2e%%\n',100*abs((p3(2)-p4(3))/p3(2)))
end