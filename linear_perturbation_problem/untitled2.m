clear
l0 = 100;
w = 0.1;
e = 2;
r = 2;
alpha = 1;
beta = 2;

x = 30:40;

H_tilde = @(x) w*x.^(alpha+1)+e*x.^(alpha) + r;
H0 = @(x) 0.1*x.^(2)+2*x.^(1) + 0.02;
func = @(x) l0*H_tilde(x)./H0(x).^3;


x1 = x(1);
xn = x(end);
Hn = H0(xn);
H1 = H0(x1);

w_co1 = l0*(xn^(1+3*beta)/Hn^3 - x1^(1+3*beta)/H1^3)/(xn-x1) ; 
e_co1 = l0*(xn^(3*beta)/Hn^3 - x1^(3*beta)/H1^3)/(xn-x1);
r_co1 = l0*(xn^(3*beta-alpha)/Hn^3 - x1^(3*beta-alpha)/H1^3)/(xn-x1);
a = w*w_co1 + e*e_co1 + r*r_co1 ;

w_co2 = l0*(xn*x1*(x1^(3*beta)/H1^3 - xn^(3*beta)/Hn^3)/(xn-x1));
e_co2 = l0*(xn*x1*(x1^(3*beta-1)/H1^3 - xn^(3*beta-1)/Hn^3)/(xn-x1));
r_co2 = l0*(xn*x1*(x1^(3*beta-alpha-1)/H1^3 - xn^(3*beta-alpha-1)/Hn^3)/(xn-x1));

b = w*w_co2 + e*e_co2 + r*r_co2 ;

interp = @(x,a,b) x.^(alpha-3*beta).*(a.*x+b);

int_dat = interp(x,a,b);

[wc,ec,rc ] = get_integral_coeff(x1,xn,H1,Hn,alpha,beta,l0);
int = w*wc+e*ec+r*rc;
disp(int)

plot(x,int_dat,'b',x, func(x),'ro',x,100./x.^4)