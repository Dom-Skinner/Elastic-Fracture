clear


load n995x846-late

nv = numel(v);
%{
% How to calculate h0_prime
h0_prime = zeros(n,numel(L));

for i = 1:numel(L)
for j = 1:n
    p1 = polyfit(lambda(end-2:end),hprime_data(nv+j,end-2:end,i).^u,1);
    h0_prime(j,i) = (l0(1)*p1(1) +p1(2)).^(1/u);
end
end

h0 = zeros(size(h0_prime));
h_coefficient_matrix = hprime_to_h_s(x,1/2,t);
for i = 1:numel(L)
h0(:,i) = h_integrate(h0_prime(:,i),x,n,t, h_coefficient_matrix,1/2);
end
%}
%plot(x(1:t-1),h0(1:t-1,1),x(1:t-1),h0(1:t-1,2),x(1:t-1),h0(1:t-1,3),x(1:t-1),h0(1:t-1,14))
dat = 1:t-1;
plot(x(dat),hprime_data(nv+dat,end,1),x(1:t-1),hprime_data(nv+dat,1,1),x(1:t-1),...
    hprime_data(nv+dat,1,13),x(1:t-1),hprime_data(nv+dat,end,13))
xlabel('x')
ylabel('h')
legend({['L = ', num2str(lambda(15))], ...
    ['L = ', num2str(lambda(16))],['L = ', num2str(lambda(17))],...
    ['L = ', num2str(lambda(18))]},'Location','NorthWest')
