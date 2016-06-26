%we analyse the forward transform properties of the scalar problem
%we pick a known function that transforms well.
%we only do this for h.

%we look at the effect of discretization on the entire function.

n = 100;
t = round(n/2);
endpoint = 20;

h = @(x) x + 1./sqrt(x);

%the functions for the kernel

K12x1int = @(x,z) real(2.*(4+(x+(-1).*z).^2).^(-2).*((-8).*x+(-1).*(4+(x+(-1).*z).^2).*z) ...
  +(1/2).*z.*log(4+(x+(-1).*z).^2)+(-1).*z.*log(x+(-1).*z));

K12sqrtx0int = @(x,z) real((1/2).*(2.*x.^(1/2).*(4+(x+(-1).*z).^2).^(-2).*(4+z.^2).^(-2).*( ...
 ...
 16.*(x+(-2).*z).*(4+z.^2)+(4+x.^2+(-2).*x.*z+z.^2).*((-3).*z.*(20+ ...
...
  z.^2)+x.*(28+z.^2)))+((sqrt(-1)*2)+(-1).*z).^(-5/2).*((-15)+(sqrt( ...
...
  -1)*(-10)).*z+2.*z.^2).*atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^( ...
...
  -1/2))+4.*z.^(-1/2).*atanh(x.^(1/2).*z.^(-1/2))+(-1).*((sqrt(-1)* ...
...
  2)+z).^(-5/2).*((-15)+(sqrt(-1)*10).*z+2.*z.^2).*atanh(x.^(1/2).*( ...
...
  (sqrt(-1)*2)+z).^(-1/2))));

%the analytic pressure
p_an = @(z) K12x1int(10^(10),z) - K12x1int(0,z) ...
    + K12sqrtx0int(10^(10),z) - K12sqrtx0int(0,z);

[x,z,soln,Pmatrix, equation_matrix, conditioning] = diagonal_matrix_test_h(n, endpoint);

%vector for input h: discretization
input = zeros(n,1);
input(1:t) = x(1:t)'.^(3/2) + ones(t,1);
input(t+1:n) = h(x(t+1:n)');

znew = tan((0.1:0.2:n-0.1)*atan(endpoint)/(n-1));

pressure = zeros(1,length(znew));

for j = 1:length(znew)
    pressure_row = diagonal_pressure_row_h(x,znew(j));
    pressure(j) = pressure_row*input;
end

plot(znew, pressure)



