clear
load n400x30-extended
K = 3*sqrt(2*pi)*KI;
s = 0.138673;
u = 4 - 6*s;
l0 = 0.5;

z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

h0=x.^(2/3)+0.5.*x.^2;
h0_z=z.^(2/3)+0.5.*z.^2;
hprime = (2/3).*x.^(-1/3)+x;
hprime(1:t-1) = (2/3).*ones(1,t-1)+x(1:t-1).^(4/3);

h_coefficient_matrix = hprime_to_h_s(x,2/3);


R2 = lubrication_integral(x,z,n,t,h_coefficient_matrix,h0,h0_z,l0,2/3);

%R = lubrication_integral_s(x,z,n,t,h_coefficient_matrix,h0,l0,0.5);
%R=R(1:n-1,n+1:2*n);

func = @(x) (x.^(2/3)+0.5*x.^2).^(-2);
int_func = @(x) integral(func,x,10000);


exact =@(x) -3.497-3/32 *(-(8*x)./(x.^(4/3)+2)-5* 2^(1/4) *log(sqrt(2) *x.^(2/3)...
    - 2^(5/4) *x.^(1/3)+2)+5*2^(1/4)*log(sqrt(2)*x.^(2/3)+2^(5/4)...
    *x.^(1/3)+2)-32./x.^(1/3)+10* 2^(1/4)*atan(1-2^(1/4)* x.^(1/3))-...
    10*2^(1/4)* atan(2^(1/4)* x.^(1/3)+1));

int = R2*hprime';
%der = (int(2:end)-int(1:end-1))./(z(2:end)'-z(1:end-1)');

%plot(x,-1./h0.^2,x(2:n-1),der,'o',)
%axis([0,12,-3.5,0.5])
for k = 1:n-1
    int_mat(k) = int_func(z(k));
end
 plot(z,int,'o',z,exact(z),z,int_mat)   
