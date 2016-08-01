%%{
n = 400;
t = 200;
xmax = 15;
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
n_increase = 8;
%}
for iterate = 1:12

scaled_K_of_c_march
K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0(iterate) = p1(2);
x_val(iterate) = xmax;

scale = x(n)/x(n-1);

for k = 1:n_increase
    x(n+k) = x(n) * scale^(k);
    z(n+k-1) = z(n-1) * scale^(k);
end

n = n+n_increase;
xmax = x(n);
end

subplot(2,1,1)
plot(x_val, l0,'o-')
xlabel('xend')
ylabel('l0')

subplot(2,1,2)
plot(x_val.^(-2), l0,'o-')
xlabel('xend^{-2}')
ylabel('l0')