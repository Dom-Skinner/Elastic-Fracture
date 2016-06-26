%we only use one of the kernels, to solve for just h.
%uses just the pressure kernel

%we look at different values of n, and of xmax, and see what happens

n = [100, 200, 400, 800];
xmax = [20, 30, 40];

soln = zeros(max(n), length(xmax), length(n));
x = zeros(max(n),length(xmax), length(n));
z = zeros(max(n),length(xmax), length(n));
condition = zeros(length(xmax),length(n));

KI = zeros(length(xmax), length(n));

for j = 1:length(xmax)
    for i = 1:length(n)
        [xtemp,ztemp,solntemp,~,~,conditiontemp] = diagonal_matrix_test_lininfty(n(i),xmax(j));
        x(1:n(i),j,i) = xtemp;
        z(1:n(i),j,i) = ztemp;
        soln(1:n(i),j,i) = solntemp;
        condition(j,i) = conditiontemp;
    end
end
     

dsoln = zeros(max(n)-1,length(xmax),length(n));
ddsoln = zeros(max(n)-2,length(xmax),length(n));

for j = 1:length(xmax)
    for i = 1:length(n)
        dsoln(1:n(i)-1,j,i) = diff(soln(1:n(i),j,i))./diff(x(1:n(i),j,i));
        ddsoln(1:n(i)-2,j,i) = diff(dsoln(1:n(i)-1,j,i))./diff(x(2:n(i),j,i));
    end
end

for j=1:length(xmax)
    for i = 1:length(n)
        a(j,i) = (soln(n(i)-1,j,i)-soln(n(i),j,i))/(x(n(i)-1,j,i)-x(n(i),j,i));
        b(j,i) = (x(n(i)-1,j,i)*soln(n(i),j,i)-x(n(i),j,i)*soln(n(i)-1,j,i))/(x(n(i)-1,j,i)-x(n(i),j,i));
    end
end

    

figure(2);
for j=1:length(xmax)
    loglog(x(3:800,j,length(n)),abs(ddsoln(1:798,j,length(n))),'.-')
    hold on;
end