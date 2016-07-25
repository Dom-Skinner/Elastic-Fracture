function [interp1, interp2] = interpolate_hprime(x,n,hprime_data,K)
% Interpolates h_0' based on values of h' for non zero K.
% Also attempts to correct for the LEFM boundary layer since this is not
% present in the h_0 solution.
%
% interp1 is the first interpolation - no LEFM correction
% interp2 is the second interpolation - LEFM corrected.
%
% interp2 is stored as usual with x^{1/2}(ax+b) for i=1,...,t-1
s = 0.138673;
u = 4 - 6*s;
interp1 = zeros(1,n);
er = zeros(1,n);
pc_er = zeros(1,n);

% First go at interpolation of hprime
for l = 1:n
    p1 = polyfit(K(end-1:end).^u , hprime_data(n+l,end-1:end),1);
    p2 = polyfit(K(end-2:end).^u , hprime_data(n+l,end-2:end),2);
    interp1(l) = p1(2);
    er(l)    = abs(p1(2)-p2(3));
    pc_er(l) = abs(p1(2)-p2(3))/p1(2);
end

fprintf('Approximate interpolation error = %.2e%%\n',100*max(pc_er(40:end)))
% Now we correct for the LEFM boundary

er = zeros(1,round(0.2*n));
for k = 1:round(0.25*n)
    p3 = polyfit([x(k:k+1), 0] , [interp1(k:k+1).*x(k:k+1).^(-1/6),0.352],1);
    p4 = polyfit([x(k:k+2), 0] , [interp1(k:k+2).*x(k:k+2).^(-1/6),0.352],2);
    er(k) = abs(p3(2)-p4(3));
end
[~,I] = min(er(:));

p3 = polyfit(x(I:I+1) , interp1(I:I+1).*x(I:I+1).^(-1/6),1);
p4 = polyfit(x(I:I+2) , interp1(I:I+2).*x(I:I+2).^(-1/6),2);
interp2 = interp1;
interp2(1:I) = (p3(1)*x(1:I)+p3(2)).*x(1:I).^(1/6);
fprintf('Approximate LEFM correction error = %.2e%%\n',100*abs((p3(2)-p4(3))/p3(2)))
end