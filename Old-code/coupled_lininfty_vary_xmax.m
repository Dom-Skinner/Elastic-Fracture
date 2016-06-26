%we only use one of the kernels, to solve for just h.
%uses just the pressure kernel

%we look at different values of n, and of xmax, and see what happens

n = [100, 140, 200, 280, 400];
xmax = [20, 40, 80];

g = zeros(max(n), length(xmax), length(n));
h = zeros(max(n), length(xmax), length(n));

x = zeros(max(n),length(xmax), length(n));
z = zeros(max(n)-1,length(xmax), length(n));
condition = zeros(length(xmax),length(n));

KI = zeros(length(xmax), length(n));
KII = zeros(length(xmax), length(n));

for j = 1:length(xmax)
    for i = 1:length(n)
        [KI(j,i), KII(j,i), xtemp,ztemp,solntemp,condition(j,i)] = linear_spline_simpleinfty(n(i),xmax(j),0,1);
        x(1:n(i),j,i) = xtemp;
        z(1:n(i)-1,j,i) = ztemp;
        h(1:n(i),j,i) = solntemp(n(i)+1:2*n(i));
        g(1:n(i),j,i) = solntemp(1:n(i));
    end
end

dg = zeros(max(n)-1,length(xmax),length(n));
ddg = zeros(max(n)-2,length(xmax),length(n));
dh = zeros(max(n)-1,length(xmax),length(n));
ddh = zeros(max(n)-2,length(xmax),length(n));

for j = 1:length(xmax)
    for i = 1:length(n)
        dg(1:n(i)-1,j,i) = diff(g(1:n(i),j,i))./diff(x(1:n(i),j,i));
        ddg(1:n(i)-2,j,i) = diff(dg(1:n(i)-1,j,i))./diff(x(2:n(i),j,i));
    end
end

for j = 1:length(xmax)
    for i = 1:length(n)
        dh(1:n(i)-1,j,i) = diff(h(1:n(i),j,i))./diff(x(1:n(i),j,i));
        ddh(1:n(i)-2,j,i) = diff(dh(1:n(i)-1,j,i))./diff(x(2:n(i),j,i));
    end
end

hold off;

figure(1);
for j=1:length(xmax)
    loglog(x(3:400,j,length(n)),abs(ddh(1:398,j,length(n))))
    legend('h```')
    hold on;
end
hold off;

figure(2);
for j=1:length(xmax)
    loglog(x(2:400,j,length(n)),abs(dg(1:399,j,length(n))))
    legend('g``')
    hold on;
end
hold off;