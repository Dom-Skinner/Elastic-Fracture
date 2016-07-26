function [kernel_matrix, interpolate_matrix] = pressure_shear_matrix_s(x,z,s)

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


    function K12_int_s = K12_s(z,ux,lx,a)
        if abs(ux - z) < 0.03 || abs(z - lx) < 0.03
            fun_non_singular = @(x1) x1.^(a).*((x1-z).^5+12.*(x1-z).^4+96.*(x1-z))./((x1-z).^2+4).^3 ;
            
            fun_remov_sing = @(theta) ( (1+theta).^a - 1)./theta ;
            if lx/z -1 > -0.95
                singular_bit = z.^(a) .*(integral(fun_remov_sing,lx/z-1,ux/z-1) + ...
                log(abs((ux-z)/(z-lx))));            
            else
                ep = 0.05;
                lp = lx/z -1;
                int_approx = integral(fun_remov_sing,-0.95,ux/z-1)+...
                    ((ep^(a+1)-lp^(a+1))/(1+a) - (ep^(a+2)-lp^(a+2))/(2+a) +...
                    (ep^(3+a)-lp^(3+a))/(3+a) -(ep^(4+a)-lp^(4+a))/(4+a))+...
                -log((1+ep)/(1+lp));
            
                singular_bit = z.^(a) .*(int_approx + ...
                log(abs((ux-z)/(z-lx))));            
                
            end
            K12_int_s = integral(fun_non_singular,lx,ux,'AbsTol',1e-5)-singular_bit;
        else
            fun_internal = @(x1) x1.^(a).*(48.*(x1-z).^2 - 64)./(x1-z)./((x1-z).^2+4).^3;
            K12_int_s = integral(fun_internal,lx,ux,'AbsTol',1e-5);
        end
    end
K12_s_fun = @K12_s;


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
%K22sqrtx1_change_this(ux,z(i))-K22sqrtx1_change_this(lx,z(i));
    function K22_int_s = K22_s(z,ux,lx,a)
        if ux > 0.1 || a > 0
            fun_internal = @(x1) -x1.^(a).*(32-24*(x1-z).^2)./((x1-z).^2+4).^3;
            K22_int_s = integral(fun_internal,lx,ux,'AbsTol',1e-5);
        else
            % rewrite an integral with integrable singularity as one
            % without, by integration by parts
            const = lx.^(a+1)./(a+1).*( 32 - 24*(lx-z).^2 ) ./ ( (lx-z).^2 + 4 ).^3 - ...
                ux.^(a+1)./(a+1).*( 32 - 24*(ux-z).^2 ) ./ ( (ux-z).^2 + 4 ).^3 ;
            fun_internal = @(x1) x1.^(a+1).*(x1-z).*( (x1-z).^2 -4) ...
                ./((x1-z).^2+4).^4;
            K22_int_s = const+(96/(a+1))*integral(fun_internal,lx,ux,'AbsTol',1e-5);
        end
    end
K22_s_fun = @K22_s;




for i=1:n-1
    %goes over the 1/sqrt(x) panels
    for j=1:t-1
        lx = x(j);
        ux = x(j+1);
        %for K11
        function_matrix(i,j) = K11sqrtx1(ux,z(i))-K11sqrtx1(lx,z(i));
        function_matrix(i,n+j) = K11sqrtx0(ux,z(i))-K11sqrtx0(lx,z(i));
        %for K12
        function_matrix(i,2*n+j) = K12_s_fun(z(i),ux,lx,s);
        %K12sqrtx1_change_this(ux,z(i))-K12sqrtx1_change_this(lx,z(i));
        function_matrix(i,2*n+n+j) = K12_s_fun(z(i),ux,lx,s-1);
        %K12sqrtx0_change_this(ux,z(i))-K12sqrtx0_change_this(lx,z(i));
        %for K21
        function_matrix(n-1+i,j) = K21sqrtx1(ux,z(i))-K21sqrtx1(lx,z(i));
        function_matrix(n-1+i,n+j) = K21sqrtx0(ux,z(i))-K21sqrtx0(lx,z(i));
        %for K22
        function_matrix(n-1+i,2*n+j) = K22_s_fun(z(i),ux,lx,s);
        function_matrix(n-1+i,2*n+n+j) = K22_s_fun(z(i),ux,lx,s-1);
        %K22sqrtx0_change_this(ux,z(i))-K22sqrtx0_change_this(lx,z(i));
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