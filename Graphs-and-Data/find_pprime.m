function pprime_data = find_pprime(lambda,x,hprime_data,n,t)
h_coefficient_matrix = hprime_to_h_l(x);

h_data = zeros(n,numel(lambda));
pprime_data = zeros(n,numel(lambda));

for l = 1:numel(lambda)
    h_data(:,l) = h_integrate(hprime_data(n+1:end,l),x,n,t,h_coefficient_matrix);
    pprime_data(:,l) = lambda(l)./(h_data(:,l).^2);
end

end