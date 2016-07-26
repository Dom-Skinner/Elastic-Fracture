clear
%load n400x40-extended
n=200;
xmax = 40;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
s=1/2;

[kernel_matrix_l, ~] = pressure_shear_matrix_l(x,z);
[kernel_matrix_s, ~] = pressure_shear_matrix_s(x,z,s);

disp(max(abs(kernel_matrix_l(:)-kernel_matrix_s(:))))
