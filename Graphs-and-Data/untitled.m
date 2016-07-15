%load n800x40.mat
num=1:99;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

%[p, h] = recover_p_h(lambda,x,z,num,hprime_data,n,t);

x1 = x(num);
z1 = z(num);

plot(z1,p(:,14),z1, -3*(4*pi*pi/243)^(1/3)*z1.^(-1/3))