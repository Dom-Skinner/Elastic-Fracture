%%{
l0_est=[];
D_est=[];
n_used = [];
%n_val = [120,160,250,310,150,200,250,300];
n_val = [160,190,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630];
flg=1;
for m = n_val
    n=m;   
    t = round(n/2);
    xmax = 15;
    x = tan((0:n-1)*atan((xmax)^(1/4))/(n-1)).^4;
    z = tan((0.5:1:n-1.5)*atan((xmax)^(1/4))/(n-1)).^4;
    scale = x(n)/x(n-1);
    
n_increase = round(0.015*m);
s = 0.138673;
u = 4 - 6*s;
l0 = [];
x_val=[];
D=[];
%}    


while x(n) < 500
    x(n+1) = x(n) * scale;
    z(n) = z(n-1) * scale;
    n = n+1;

end
xmax = x(n);

upper_lim = 1000;
while xmax < upper_lim
    if x(n) * scale^(n_increase) < upper_lim
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

x_m1 = x_val.^(-1);

p1 = polyfit(x_m1,D,1);
p2 = polyfit(x_m1,l0,1);
D_est(end+1) = p1(2);
l0_est(end+1) = mean(l0);%p2(2);
n_used(end+1) = n;

end

subplot(2,2,1)
plot(n_used.^(-2),D_est,'o-')
xlabel('n^{-2}')
ylabel('D')

subplot(2,2,2)
plot(n_used.^(-0.67),l0_est,'o-')
xlabel('n^{-0.5}')
ylabel('l0')

subplot(2,2,3)
plot(n_used.^(-1),D_est,'o-')
xlabel('n^{-1}')
ylabel('D')

subplot(2,2,4)
plot(n_used.^(-1),l0_est,'o-')
xlabel('n^{-1}')
ylabel('l0')
