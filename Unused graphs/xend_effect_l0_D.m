l0_est=[];
D_est=[];
n_used = [];
n_val = [200,300,400,500]; 

for m = n_val
    n=m;
   
    t = round(n/2);

    xmax = 15;
    x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
    z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
    scale = x(n)/x(n-1);
    
n_increase = round(0.013*m);
s = 0.138673;
u = 4 - 6*s;
l0 = [];
x_val=[];
D=[];


while x(n) < 400
    x(n+1) = x(n) * scale;
    z(n) = z(n-1) * scale;
    n = n+1;

end
xmax = x(n);

upper_lim = 1000;
while xmax < upper_lim
    lambda = 0.0584:0.0004:0.0588;
    
    scaled_K_of_c_march
    K = 3*sqrt(2*pi)*KI;
    p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
    % The error in this interpolation has been tested to be around 10^-6
    % and so, should not pose a problem
    %p2 = polyfit(K(end-2:end).^u,lambda(end-2:end),2);
    %disp(abs(p1(2)-p2(3)))
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
%p2 = polyfit(x_m1,D,2); %
%disp(abs(p1(2)-p2(3))) %
p2 = polyfit(x_m1,l0,1);
%p1 = polyfit(x_m1,l0,2); % Despite our poor understanding and real
%disp(abs(p1(3)-p2(2))) % inability to extrapolate here, the error is maybe
% 10^-6, certainly an order of magnitude less than the error due to the tip
% resolution.
D_est(end+1) = p1(2);
l0_est(end+1) = mean(l0);
n_used(end+1) = n;

hold on
subplot(2,1,1)
plot(x_m1, D,'o-',[x_m1 0],[x_m1 0]*p1(1) +p1(2),':')  


xlabel('xend^{-1}')
ylabel('D')

hold on
subplot(2,1,2)
plot(x_m1, l0,'o-',[x_m1 0],[x_m1 0]*p2(1) +p2(2),':')
xlabel('xend^{-1}')
ylabel('l0')
end