clear 
load n800x40.mat

num = 1:99;
s = 0.138673;
u = 4 - 6*s;
y1 = x(num)';
lend = numel(lambda);

h = zeros(numel(num),lend);
h_coefficient_matrix = hprime_to_h(x);

for l = 1:14
    hprime = hprime_data(:,l);
    h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);
    h_temp = zeros(1,t);
    for r = 1:t-1
        h_temp(r) = h_coeffs(3*r -2)*x(r)^(3/2) + ...
            h_coeffs(3*r-1)*x(r)^(1/2) + h_coeffs(3*r);
    end
    h(:,l) = lambda(l).^(-1/3).*h_temp(num);
end

z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
[~, h] = recover_p_h(lambda,x,z,num,hprime_data,n,t);

K = 3*sqrt(2*pi)*lambda.^(-1/3).*KI;
Kap = 3*sqrt(2*pi)*KI;


% interpolate to find h_0
h0 = (K(lend)^u*h(:,lend-1)-K(lend-1)^u*h(:,lend))/(K(lend)^u-K(lend-1)^u); 
% interpolate to find lambda_0. Note we use KI here
l0 =  (Kap(lend)^u*lambda(lend-1)-Kap(lend-1)^u*lambda(lend))/...
    (Kap(lend)^u-Kap(lend-1)^u);
ch = zeros(1,numel(num));
h02 = zeros(1,numel(num));
for l = num
    p = polyfit(K(8:14).^u,h(l,8:14),1);
    h02(l)=p(2);
    ch(l) = p(1);
end
%plot(y1.^s, ch)
p = polyfit(y1.^s,ch',1);
plot(y1.^s, ch,'*-',y1.^s, p(1).*y1.^s+p(2))

return

% Estimate the value of C, by finding the value which minimises the error
er=10;
for w = 0.002:0.0001:0.004
    er1 = 0;
    for l = 9:11
            er1 = er1+ norm( h0 + y1.^(s)*K(l)^u*w - h(:,l));
    end
    if er > er1
        er = er1;
        c = w;
    end
end
fprintf('C = %2.1e\n', c);
fprintf('lambda_0 = %2.2e\n', lambda(14));
