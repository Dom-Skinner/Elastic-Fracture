%this is a check to see if our solution works!!! by confirming it for the
%earlier work

Mload = 1;
Pload = 0;

n = 400;
t = round(n/2);
xmax = 30;

%the data points
%x = tan((0:n-1)*atan(xmax)/(n-1));
x(1:t) = tan((0:t-1)*(pi/4)/(t-1));
x(t+1:n) = (1:n-t)*(xmax-1)/(n-t) + ones(1,n-t);

%z = tan((0.5:1:n-1.5)*atan(xmax)/(n-1));
z(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1));
z(t:n-1) = (0.5:n-t-0.5)*(xmax-1)/(n-t) + ones(1,n-t);

[kernel_matrix, interpolate_matrix] = pressure_shear_matrix(x,z);

gprime_infty = interpolate_matrix(2*n,:);
hprime_infty = interpolate_matrix(3*n,:);

equation_matrix = zeros(2*n,2*n);
equation_matrix(1:2*n-2,:) = kernel_matrix;
equation_matrix(2*n-1,:) = gprime_infty;
equation_matrix(2*n,:) = hprime_infty;

forcing = zeros(2*n,1);
forcing(2*n-1) = (-Pload+6*Mload);
forcing(2*n) = (12*Mload);

%normalizes
forcing = (2*pi)^(3/2)/(8*pi)*forcing;

solution = equation_matrix\forcing;

KI = solution(n+1);
KII = solution(1);