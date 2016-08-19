clear 
load xend_march
% The code that created the data in xend_march is commented out below.

%{
l0_est=[];
D_est=[];
n_used=[];
   
n = 400;
t = round(n/2);

xmax = 15;
x = tan((0:n-1)*atan((xmax)^(1/3))/(n-1)).^3;
z = tan((0.5:1:n-1.5)*atan((xmax)^(1/3))/(n-1)).^3;
scale = x(n)/x(n-1);
    
n_increase = round(0.015*n);
s = 0.138673;
u = 4 - 6*s;
l0 = [];
x_val=[];
D=[];
    


while x(n) < 10
    x(n+1) = x(n) * scale;
    z(n) = z(n-1) * scale;
    n = n+1;

end
xmax = x(n);

upper_lim = 600;
while xmax < upper_lim
    if x(n) * scale^(n_increase) <1950
        lambda = 0.0584:0.0004:0.0588;
    else
        %lambda = 0.056:0.0004:0.0588;
    end
    scaled_K_of_c_march
    K = 3*sqrt(2*pi)*KI;
    p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
    l0(end+1) = p1(2);
    D(end+1) = p1(1);
    x_val(end+1) = xmax;

    if x(n) * scale^(n_increase) < upper_lim
        for k = 1:n_increase
            x(n+k) = x(n) * scale^(k);
            z(n+k-1) = z(n-1) * scale^(k);
        end

        n = n+n_increase;
        xmax = x(n);
    
    else
        break
    end
end

%}

hold on
subplot(1,2,1)
plot(x_val.^(-1),D,'o-')
xlabel('x^{-1}')
ylabel('D')

hold on
subplot(1,2,2)
plot(x_val.^(-1),l0,'o-')
xlabel('x^{-1}')
ylabel('l0')
%axis([0,0.0016,0.05908,0.05918])
%}