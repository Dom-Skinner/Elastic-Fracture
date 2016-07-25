clear

w = 0.1;
e = 2;
r = 2;
alpha = 0.14;
beta = 2/3;

x = 0.003:0.001:0.05;

H_tilde = @(x) w*x.^(alpha+1)+e*x.^(alpha) + r;
H0 = @(x) 0.1*x.^(beta+1)+2*x.^(beta) + 0.02;


xmeas = x(1:5:end);
for k = 1:numel(xmeas)-1
x1 = xmeas(k);
xn = xmeas(k+1);
Hn = H0(xn);
H1 = H0(x1);

a(k) = w*(xn^(1+2*beta)/Hn^2 - x1^(1+2*beta)/H1^2)/(xn-x1) + ...
    e*(xn^(2*beta)/Hn^2 - x1^(2*beta)/H1^2)/(xn-x1) + ...
    r*(xn^(2*beta-alpha)/Hn^2 - x1^(2*beta-alpha)/H1^2)/(xn-x1) ;

b(k) = w*xn*x1*(x1^(2*beta)/H1^2 - xn^(2*beta)/Hn^2)/(xn-x1) + ...
    e*xn*x1*(x1^(2*beta-1)/H1^2 - xn^(2*beta-1)/Hn^2)/(xn-x1) + ...
    r*xn*x1*(x1^(2*beta-alpha-1)/H1^2 - xn^(2*beta-alpha-1)/Hn^2)/(xn-x1) ;

al(k) = (H_tilde(xn)/H0(xn)^2-H_tilde(x1)/H0(x1)^2)/(xn-x1);
bl(k) = H_tilde(xn)/H0(xn)^2-al(k)*xn;
end

interp = @(x,a,b) x.^(alpha-2*beta).*(a.*x+b);
interpl = @(x,a,b) (a.*x+b);

hold on
int = 0;
for k = 1:numel(xmeas)-1
    int_dat = interp(x(5*k-4:5*k+1),a(k),b(k));
    plot(x(5*k-4:5*k+1),int_dat,'b');
    int_datl = interpl(x(5*k-4:5*k+1),al(k),bl(k));
    plot(x(5*k-4:5*k+1),int_datl,'k');
    int = int + [w,e,r].*get_integral_coeff(x(5*k-4),x(5*k+1),int_dat(1),int_dat(6),alpha,beta,1);
end

plot(x, H_tilde(x)./H0(x).^2,'ro')