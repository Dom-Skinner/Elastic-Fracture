%%{
l0_est=[];
D_est=[];
n_val = [200,250,300,350,400,500,600,750];
for m = n_val
    n=m;
t = round(n/2);
xmax = 15;
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
n_increase = round(0.02*m);
s = 0.138673;
u = 4 - 6*s;
l0 = [];
x_val=[];
D=[];
%}    
scale = x(n)/x(n-1);

while x(n) < 400
    x(n+1) = x(n) * scale;
    z(n) = z(n-1) * scale;
    n = n+1;

end
xmax = x(n);

while xmax < 850
    if x(n) * scale^(n_increase) <850
        lambda = 0.0584:0.0004:0.0588;
    else
        lambda = 0.056:0.0004:0.0588;
    end
    scaled_K_of_c_march
    K = 3*sqrt(2*pi)*KI;
    p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
    l0(end+1) = p1(2);
    D(end+1) = p1(1);
    x_val(end+1) = xmax;

    if x(n) * scale^(n_increase) <850
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

x_m2 = x_val.^(-2);

p1 = polyfit(x_m2,D,1);
p2 = polyfit(x_m2,l0,1);
D_est(end+1) = p1(2);
l0_est(end+1) = p2(2);

%%% Want to save:
file = ['n' , num2str(n), 'x', num2str(round(x(n)))];
D = D_est(end);
l0 = l0_est(end);
save(file,'D','l0','n','hprime_data','lambda','KI','t','x','z')
%{
hold on
subplot(2,1,1)
plot(x_m2, D,'o-',[x_m2 0],[x_m2 0]*p1(1) +p1(2),':')
xlabel('xend^{-2}')
ylabel('D')

hold on
subplot(2,1,2)
plot(x_m2, l0,'o-',[x_m2 0],[x_m2 0]*p2(1) +p2(2),':')
xlabel('xend^{-2}')
ylabel('l0')
%}
end
%legend({'200','300','400','500'})
%w = waitforbuttonpress;
hold off
subplot(2,1,1)
plot(n_val.^(-1),D_est,'o-')
xlabel('n^{-1}')
ylabel('D')
axis([0,0.005,-0.0083,-0.00806])

hold off
subplot(2,1,2)
plot(n_val.^(-1),l0_est,'o-')
xlabel('n^{-1}')
ylabel('l0')
axis([0,0.005,0.0591,0.0594])