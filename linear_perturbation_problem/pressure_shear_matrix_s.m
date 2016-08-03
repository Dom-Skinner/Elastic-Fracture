function [kernel_matrix, interpolate_matrix] = pressure_shear_matrix_s(x,z,s,t)

n = length(x);
infinity = 10^10;

interpolate_matrix = linear_simpleinfty_interpolate_s(x,s);
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
        fun0 = @(x1) x1.^(a+1).*((x1-z).^5+12.*(x1-z).^3+96.*(x1-z)) ... 
            ./((x1-z).^2+4).^3/(a+1) ;
        fun1 = @(x1) x1.^(a).*((x1-z).^5+12.*(x1-z).^3+96.*(x1-z))./((x1-z).^2+4).^3 ;
        fun2 = @(x1) x1.^(a+1)./(a+1)./( (x1-z).^2 + 4).^4 .*( (x1-z).^6 + ...
            16.*(x1-z).^4 +336.*(x1-z).^2 - 384 );
        fun3 = @(theta) ((1+theta).^a-1)./theta;
        fun4 = @(u) u.^(a+1)./(u-1).^2 /(a+1);
        fun5 = @(u) u.^(a+1)./(u-1) /(a+1);
        
        if a < 0 && lx < 0.5
            non_singular_bit = fun0(ux) -fun0(lx) + ...
                integral(fun2, lx,ux,'AbsTol',1e-5);
        else
            non_singular_bit = integral(fun1,lx,ux,'AbsTol',1e-5);
        end
        
        if a >= 0 ||  lx/z > 0.5
            singular_bit = -z^a *( integral(fun3, lx/z -1, ux/z -1, ...
                'AbsTol',1e-5) + log(abs((ux-z)/(z-lx))));
        elseif ux/z < 0.5 
            singular_bit = z^a *( fun5(lx/z)-fun5(ux/z)  - ...
                integral(fun4,lx/z,ux/z,'AbsTol',1e-5) );
        else
               singular_bit = z^a *( (lx/z)^(a+1)/(a+1)/((lx/z)-1) + ...
                (2)^(-a)/(a+1) - ...
                integral(fun4,lx/z,0.5,'AbsTol',1e-5) - ...
                integral(fun3,-0.5,ux/z-1,'AbsTol',1e-5)-log(2*abs(ux/z-1))) ;
        
        end
        K12_int_s = singular_bit + non_singular_bit;
        return 
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
        if  ux > 0.1 || a > 0
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
        x_t = x(1:t);
    %
    %atan_store_1 = atan(x_t.^(1/2).*((sqrt(-1)*2)+(-1)*z(i)).^(-1/2));
    %atan_store_2 = atan(x_t.^(1/2).*((sqrt(-1)*(-2))+(-1)*z(i)).^(-1/2));
    %
   % atanh_store_1 = atanh(x_t.^(1/2).*((sqrt(-1)*2)+z(i)).^(-1/2));
   % atanh_store_2 = atanh(x_t.^(1/2).*z(i).^(-1/2));
    % for K11
    K11sqrtx1_store = K11sqrtx1(x_t,z(i));%,atan_store_1,atan_store_2); 
    K11sqrtx0_store = K11sqrtx0(x_t,z(i));%,atan_store_1,atan_store_2); 
    % for K21
    K21sqrtx1_store = K21sqrtx1(x_t,z(i));%,atan_store_1,atanh_store_1,atanh_store_2);
    K21sqrtx0_store = K21sqrtx0(x_t,z(i));%,atan_store_1,atanh_store_1,atanh_store_2);
   
    
    j=1:t-1;
    %for K11
    function_matrix(i,j)         = K11sqrtx1_store(j+1)-K11sqrtx1_store(j);
    function_matrix(i,n+j)       = K11sqrtx0_store(j+1)-K11sqrtx0_store(j);
    %for K21
    function_matrix(n-1+i,j)     = K21sqrtx1_store(j+1)-K21sqrtx1_store(j);
    function_matrix(n-1+i,n+j)   = K21sqrtx0_store(j+1)-K21sqrtx0_store(j);
    
    
    for j=1:t-1
        lx = x(j);
        ux = x(j+1);
        %for K12
        function_matrix(i,2*n+j) = K12_s_fun(z(i),ux,lx,s);
        function_matrix(i,2*n+n+j) = K12_s_fun(z(i),ux,lx,s-1);
     
        %for K22
        function_matrix(n-1+i,2*n+j) = K22_s_fun(z(i),ux,lx,s);
        function_matrix(n-1+i,2*n+n+j) = K22_s_fun(z(i),ux,lx,s-1);
    end
    %goes over the linear panels including infinity
    x_tmp = [x(t:n) infinity];
    K11x1_store = K11x1(x_tmp,z(i));
    K11x0_store = K11x0(x_tmp,z(i));
    K12x1_store = K12x1(x_tmp,z(i));
    K12x0_store = K12x0(x_tmp,z(i));
    K21x1_store = K21x1(x_tmp,z(i));
    K21x0_store = K21x0(x_tmp,z(i));
    K22x1_store = K22x1(x_tmp,z(i));
    K22x0_store = K22x0(x_tmp,z(i));
    %
    j=1:n-t+1;
    %
    function_matrix(i,j+t-1)           = K11x1_store(j+1)-K11x1_store(j);
    function_matrix(i,n+j+t-1)         = K11x0_store(j+1)-K11x0_store(j);
    %for K12
    function_matrix(i,2*n+j+t-1)       = K12x1_store(j+1)-K12x1_store(j);
    function_matrix(i,2*n+n+j+t-1)     = K12x0_store(j+1)-K12x0_store(j);
    %for K21
    function_matrix(n-1+i,j+t-1)       = K21x1_store(j+1)-K21x1_store(j);
    function_matrix(n-1+i,n+j+t-1)     = K21x0_store(j+1)-K21x0_store(j);
    %for K22
    function_matrix(n-1+i,2*n+j+t-1)   = K22x1_store(j+1)-K22x1_store(j);
    function_matrix(n-1+i,2*n+n+j+t-1) = K22x0_store(j+1)-K22x0_store(j);
    disp(i)
end        

kernel_matrix = function_matrix*interpolate_matrix;

return
end