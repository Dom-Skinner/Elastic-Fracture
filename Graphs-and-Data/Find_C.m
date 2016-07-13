num = 1:99;
s = 0.138673;
y1 = x(num)';

h = zeros(numel(num),numel(lambda));
h_coefficient_matrix = hprime_to_h(x);

for l = 1:14
hprime = hprime_data(:,l);
h_coeffs = h_coefficient_matrix*hprime(n+1:2*n);
h_temp = zeros(1,t);
for r = 1:t-1
    h_temp(r) = h_coeffs(3*r -2)*x(r)^(3/2) + h_coeffs(3*r-1)*x(r)^(1/2) + ...
        h_coeffs(3*r);
end
h(:,l) = lambda(l).^(-1/3).*h_temp(num);
end



K = 3*sqrt(2*pi)*lambda.^(-1/3).*KI;
%disp(K)

for m = 1:-1
%    x1 = lambda(m)^(-1/3).*hprime_data(n+num,m);
%    x2 = lambda(14)^(-1/3).*hprime_data(n+num,14);
%    diff = (h(:,m)-h(:,14))./(K(m)^(4-6*s)-K(14)^(4-6*s));
%    plot(y1, diff.*y1.^(-s));
    plot(y1, h(:,m), '*-', y1, (2*K(m)/sqrt(18*pi))*y1.^(1/2), ...
        y1, (243/(4*pi*pi))^(1/6) *y1.^(2/3));
    s = '$K_I$ =';
    s1 = strcat(s, num2str(K(m)));
    title(s1,'Interpreter','latex','fontsize',20')
    legend({'$h$','$Kx^{1/2}$','$A_0 x^{1/3}$'},'Interpreter','latex','fontsize',20')
    axis([0,0.02,0,0.2])
    
    hold off
    %hold on
    pause(3)
end

% Estimate the value of C, by finding the value which minimises the error
er=10;
for w = 0.002:0.0001:0.004
    er1 = 0;
    for l = 9:11
            er1 = er1+ norm( h(:,14) + y1.^(s)*K(l)^(4-6*s)*w - h(:,l));
    end
    if er > er1
        er = er1;
        c = w;
    end
end
fprintf('C = %2.1e\n', c);
fprintf('lambda_0 = %2.2e\n', lambda(14));

% Same trick with estimating $\lambda_1$

er=10;
for w = -3:0.001:-1
    er1 = 0;
    for l = 9:11
            er1 = er1+ norm( lambda(l)-lambda(14)-c*lambda(14)^(1-2*s)*K(l)^(4-6*s)*w);
    end
    if er > er1
        er = er1;
        l1 = w;
    end
end
fprintf('lambda_1 = %2.2e\n', l1);

plot(K,lambda,K(9:14), lambda(14)+c*l1*K(9:14).^(4-6*s)*lambda(14)^(1-2*s))

%for m = 9
%    plot(y1,h(:,m),'*-',y1,h(:,14), y1,h(:,14) + y1.^(s)*K(m)^(4-6*s)*0.0028);
%end
%plot(y1, (243/(4*pi*pi))^(1/6) *y1.^(2/3));
%plot(y1, (2*K(6)/sqrt(18*pi))*y1.^(1/2));