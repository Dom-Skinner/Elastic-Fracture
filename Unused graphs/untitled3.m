
x_vals = [20,25,30,40];

%n=150;

for k = 1:numel(x_vals)
    xmax = x_vals(k);
    scaled_K_of_c_march
    s = 0.138673;
    u = 4 - 6*s;
    K = 3*sqrt(2*pi)*KI;
    p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
    l0(k) = p1(2);
end
hold on
plot(x_vals.^(-2),l0,'o-')