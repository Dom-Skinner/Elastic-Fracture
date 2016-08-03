clear
%load n273x778
n = 40;
t=20;
s = 0.1387;
xmax=20;
x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

K1 = pressure_shear_matrix_s(x,z,s,t);
K2 = pressure_shear_matrix_s_vectorised(x,z,s,t);

disp(max(abs(K1(:)-K2(:))))