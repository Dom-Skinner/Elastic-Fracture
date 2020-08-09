function pprime_data = find_pprime(lambda,x,hprime_data,n,t)
% Finds the derivative of the pressure by first integrating hprime, and
% then using the lubrication formula.
h_coefficient_matrix = hprime_to_h_s(x,0.5,t);

h_data = zeros(n,numel(lambda));
pprime_data = zeros(n,numel(lambda));

for l = 1:numel(lambda)
    h_data(:,l) = h_integrate(hprime_data(n+1:end,l),x,n,t,h_coefficient_matrix,0.5);
    pprime_data(:,l) = lambda(l)./(h_data(:,l).^2);
end

end