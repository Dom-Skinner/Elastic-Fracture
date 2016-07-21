function p0_prime = find_p0_prime(n,pprime_data,K)
% Interpolates p_0' based on values of p' for non zero K.
% Also attempts to correct for the LEFM boundary layer since this is not
% present in the p_0 solution.
%
%
s = 0.138673;
u = 4 - 6*s;
p0_prime = zeros(1,n);
er = zeros(1,n);
pc_er = zeros(1,n);

% First go at interpolation of pprime
for l = 1:n
    p1 = polyfit(K(end-1:end).^u , pprime_data(l,end-1:end),1);
    p2 = polyfit(K(end-2:end).^u , pprime_data(l,end-2:end),2);
    p0_prime(l) = p1(2);
    er(l)    = abs(p1(2)-p2(3));
    pc_er(l) = abs(p1(2)-p2(3))/p1(2);
end

fprintf('Approximate interpolation error = %.2e%%\n',100*max(pc_er(40:end)))
end