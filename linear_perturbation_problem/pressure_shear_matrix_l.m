function [kernel_matrix, interpolate_matrix] = pressure_shear_matrix_l(x,z)

n = length(x);
t = round(n/2);
infinity = 10^10;

interpolate_matrix = linear_simpleinfty_interpolate_s(x,0.5);
function_matrix = zeros(2*(n-1),4*n);


%FINDS ALL THE KERNEL FUNCTIONS

K11x1 = @(x,z) real(4.*(4+3.*x.^2+(-4).*x.*z+z.^2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2));
    
K11x0 = @(x,z) real(8.*(4+(x+(-1).*z).^2).^(-2).*(x+(-1).*z));


%finds the kernel K12

K12x1 = @(x,z) real(2.*(4+(x+(-1).*z).^2).^(-2).*((-8).*x+(-1).*(4+(x+(-1).*z).^2).*z) ...
  +(1/2).*z.*log(4+(x+(-1).*z).^2)+(-1).*z.*log(x+(-1).*z));
  
  
K12x0 = @(x,z) real((-2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*(12+x.^2+(-2).*x.*z+z.^2)+( ...
  1/2).*log(4+(x+(-1).*z).^2)+(-1).*log(x+(-1).*z));


%finds the kernel K21

K21x1 = @(x,z) -real(16.*x.*(4+(x+(-1).*z).^2).^(-2)+(-8).*x.*(4+(x+(-1).*z).^2).^(-1)+ ...
  2.*(4+(x+(-1).*z).^2).^(-1).*z+4.*atan((1/2).*(x+(-1).*z))+(-1/2) ...
  .*z.*log(4+(x+(-1).*z).^2)+z.*log(x+(-1).*z));
  
    
K21x0 = @(x,z) -real(2.*(8+(-3).*(4+(x+(-1).*z).^2)).*(4+(x+(-1).*z).^2).^(-2)+(-1/2).* ...
  log(4+(x+(-1).*z).^2)+log(x+(-1).*z));


%finds the kernel K22

K22x1 = @(x,z) -real(4.*(4+3.*x.^2+(-4).*x.*z+z.^2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2));

    
K22x0 = @(x,z) -real(8.*(4+(x+(-1).*z).^2).^(-2).*(x+(-1).*z));

K11sqrtx1 = @(x,z) real((-2).*x.^(1/2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*( ...
...
  x.^3+(-4).*x.^2.*z+x.*((-12)+z.^2)+2.*z.*(4+z.^2))+(sqrt(-1)*(1/2) ...
...
  ).*((sqrt(-1)*(-2))+(-1).*z).^(-3/2).*atan(x.^(1/2).*((sqrt(-1)*( ...
...
  -2))+(-1).*z).^(-1/2))+(sqrt(-1)*(-1/2)).*((sqrt(-1)*2)+(-1).*z) ...
...
  .^(-3/2).*atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2)));
  
    
K11sqrtx0 = @(x,z) real(2.*x.^(1/2).*(4+z.^2).^(-2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*(6.* ...
...
  x.^3.*z+x.^2.*(4+(-23).*z.^2)+32.*x.*(z+z.^3)+(-5).*((-16)+8.* ...
 ...
 z.^2+3.*z.^4))+(sqrt(-1)*(3/2)).*((sqrt(-1)*(-2))+(-1).*z).^(-5/2) ...
...
  .*atan(x.^(1/2).*((sqrt(-1)*(-2))+(-1).*z).^(-1/2))+(sqrt(-1)*( ...
...
  -3/2)).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*atan(x.^(1/2).*((sqrt(-1)* ...
...
  2)+(-1).*z).^(-1/2)));

%finds the kernel K12

K12sqrtx1 = @(x,z) real((1/2).*(2.*x.^(1/2).*(4+(x+(-1).*z).^2).^(-2).*((-16)+((-4)+x.*z+( ...
...
  -3).*z.^2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2))+(-1).*(( ...
 ...
 sqrt(-1)*2)+(-1).*z).^(-3/2).*((-3)+(sqrt(-1)*(-6)).*z+2.*z.^2).* ...
...
  atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+4.*z.^(1/2).*atanh( ...
...
  x.^(1/2).*z.^(-1/2))+(-1).*((sqrt(-1)*2)+z).^(-3/2).*((-3)+(sqrt( ...
...
  -1)*6).*z+2.*z.^2).*atanh(x.^(1/2).*((sqrt(-1)*2)+z).^(-1/2))));
  
  
K12sqrtx0 = @(x,z) real((1/2).*(2.*x.^(1/2).*(4+(x+(-1).*z).^2).^(-2).*(4+z.^2).^(-2).*( ...
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



%finds the kernel K21

K21sqrtx1 = @(x,z) -real((1/2).*(32.*x.^(1/2).*((4+(x+(-1).*z).^2).^(-2)+(-1/16).*(4+z.^2) ...
...
  .^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-1).*(28+x.*z+5.*z.^2))+((sqrt( ...
...
  -1)*2)+(-1).*z).^(-3/2).*((-11)+(sqrt(-1)*(-10)).*z+2.*z.^2).* ...
 ...
 atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(-4).*z.^(1/2).* ...
 ...
 atanh(x.^(1/2).*z.^(-1/2))+((sqrt(-1)*2)+z).^(-3/2).*((-11)+(sqrt( ...
...
  -1)*10).*z+2.*z.^2).*atanh(x.^(1/2).*((sqrt(-1)*2)+z).^(-1/2))));
    
  
K21sqrtx0 = @(x,z) -real((1/2).*(2.*x.^(1/2).*(4+(x+(-1).*z).^2).^(-2).*(4+z.^2).^(-2).*(( ...
...
  -16).*(x+(-2).*z).*(4+z.^2)+(4+x.^2+(-2).*x.*z+z.^2).*(4.*x+(-4).* ...
...
  z+7.*x.*z.^2+(-13).*z.^3))+(-1).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*( ...
...
  (-7)+(sqrt(-1)*(-6)).*z+2.*z.^2).*atan(x.^(1/2).*((sqrt(-1)*2)+( ...
...
  -1).*z).^(-1/2))+(-4).*z.^(-1/2).*atanh(x.^(1/2).*z.^(-1/2))+(( ...
...
  sqrt(-1)*2)+z).^(-5/2).*((-7)+(sqrt(-1)*6).*z+2.*z.^2).*atanh(x.^( ...
...
  1/2).*((sqrt(-1)*2)+z).^(-1/2))));


%finds the kernel K22

K22sqrtx1 = @(x,z) -real((-2).*x.^(1/2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*( ...
...
  x.^3+(-4).*x.^2.*z+x.*((-12)+z.^2)+2.*z.*(4+z.^2))+(sqrt(-1)*(1/2) ...
...
  ).*((sqrt(-1)*(-2))+(-1).*z).^(-3/2).*atan(x.^(1/2).*((sqrt(-1)*( ...
...
  -2))+(-1).*z).^(-1/2))+(sqrt(-1)*(-1/2)).*((sqrt(-1)*2)+(-1).*z) ...
...
  .^(-3/2).*atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2)));
  
    
K22sqrtx0 = @(x,z) -real(2.*x.^(1/2).*(4+z.^2).^(-2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*(6.* ...
...
  x.^3.*z+x.^2.*(4+(-23).*z.^2)+32.*x.*(z+z.^3)+(-5).*((-16)+8.* ...
 ...
 z.^2+3.*z.^4))+(sqrt(-1)*(3/2)).*((sqrt(-1)*(-2))+(-1).*z).^(-5/2) ...
...
  .*atan(x.^(1/2).*((sqrt(-1)*(-2))+(-1).*z).^(-1/2))+(sqrt(-1)*( ...
...
  -3/2)).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*atan(x.^(1/2).*((sqrt(-1)* ...
...
  2)+(-1).*z).^(-1/2)));




for i=1:n-1
    %goes over the 1/sqrt(x) panels
    for j=1:t-1
        lx = x(j);
        ux = x(j+1);
        %for K11
        function_matrix(i,j) = K11sqrtx1(ux,z(i))-K11sqrtx1(lx,z(i));
        function_matrix(i,n+j) = K11sqrtx0(ux,z(i))-K11sqrtx0(lx,z(i));
        %for K12
        function_matrix(i,2*n+j) = K12sqrtx1(ux,z(i))-K12sqrtx1(lx,z(i));
        function_matrix(i,2*n+n+j) = K12sqrtx0(ux,z(i))-K12sqrtx0(lx,z(i));
        %for K21
        function_matrix(n-1+i,j) = K21sqrtx1(ux,z(i))-K21sqrtx1(lx,z(i));
        function_matrix(n-1+i,n+j) = K21sqrtx0(ux,z(i))-K21sqrtx0(lx,z(i));
        %for K22
        function_matrix(n-1+i,2*n+j) = K22sqrtx1(ux,z(i))-K22sqrtx1(lx,z(i));
        function_matrix(n-1+i,2*n+n+j) = K22sqrtx0(ux,z(i))-K22sqrtx0(lx,z(i));
    end
    %goes over the linear panels including infinity
    for j=t:n
        if j==n
            lx = x(n);
            ux = infinity;
        else
            lx = x(j);
            ux = x(j+1);
        end
        function_matrix(i,j) = K11x1(ux,z(i))-K11x1(lx,z(i));
        function_matrix(i,n+j) = K11x0(ux,z(i))-K11x0(lx,z(i));
        %for K12
        function_matrix(i,2*n+j) = K12x1(ux,z(i))-K12x1(lx,z(i));
        function_matrix(i,2*n+n+j) = K12x0(ux,z(i))-K12x0(lx,z(i));
        %for K21
        function_matrix(n-1+i,j) = K21x1(ux,z(i))-K21x1(lx,z(i));
        function_matrix(n-1+i,n+j) = K21x0(ux,z(i))-K21x0(lx,z(i));
        %for K22
        function_matrix(n-1+i,2*n+j) = K22x1(ux,z(i))-K22x1(lx,z(i));
        function_matrix(n-1+i,2*n+n+j) = K22x0(ux,z(i))-K22x0(lx,z(i));        
    end    
end        

kernel_matrix = function_matrix*interpolate_matrix;

return
end
    