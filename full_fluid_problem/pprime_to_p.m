%takes in coefficients for h, stores as a 3*n vector, and gives out the
%values of p at the z values

function pressure = pprime_to_p(x, z, h_coeffs,t)

n = length(x);
infinity = 10^(10);

a = h_coeffs(1:3:3*n-2);
b = h_coeffs(2:3:3*n-1);
c = h_coeffs(3:3:3*n);


h_quad_int = @(x,a,b,c) (-1).*(b.^2+(-4).*a.*c).^(-1).*((b+2.*a.*x).* ...
  (c+x.*(b+a.*x)).^(-1)+4.*a.*((-1).*b.^2+4.*a.*c).^(-1/2).* ...
  atan(((-1).*b.^2+4.*a.*c).^(-1/2).*(b+2.*a.*x)));

%does the h function

presum_sqrt_int = @(x,a,b,c,xi) ((4.*b.^2.* ...
  log(x.^(1/2)+(-1).*xi)+9.*a.*c.*log(x.^(1/2)+(-1).*xi).* ...
  xi).*(b+3.*a.*xi.^2).^(-1));

roots_sqrt_int = @(a,b,c) roots([a,0,b,c]);

h_sqrt_int = @(x,a,b,c) 2.*(4.*b.^3+27.*a.*c.^2).^(-1).*(c+x.^(1/2).* ...
    (b+a.*x)).^(-1).*(6.*b.*c+2.*b.^2.*x.^(1/2)+9.*a.*c.*x+(c+ ...
    x.^(1/2).*(b+a.*x)).*sum(presum_sqrt_int(x,a,b,c, ...
    roots_sqrt_int(a,b,c))));

%now we actually calculate the pressure at z(j)
pressure = zeros(n-1,1);

for j=1:n-1
    for i=j:n
        if i==j
            if j <= t-1
                pressure(j) = pressure(j) + h_sqrt_int(x(j+1),a(j), ... 
                    b(j),c(j)) - h_sqrt_int(z(j),a(j),b(j),c(j));
            else
                pressure(j) = pressure(j) + h_quad_int(x(j+1),a(j), ...
                    b(j),c(j)) - h_quad_int(z(j),a(j),b(j),c(j));
            end
        elseif i~=n
            if i<=t-1
                pressure(j) = pressure(j) + h_sqrt_int(x(i+1),a(i), ...
                    b(i),c(i)) - h_sqrt_int(x(i),a(i),b(i),c(i));
            else
                pressure(j) = pressure(j) + h_quad_int(x(i+1),a(i), ... 
                    b(i),c(i)) - h_quad_int(x(i),a(i),b(i),c(i));
            end                
        else
            pressure(j) = pressure(j) + h_quad_int(infinity,a(n), ...
                b(n),c(n)) - h_quad_int(x(n),a(n),b(n),c(n));
        end
    end
end

%oops we computed \int_z^\infty, not \int_infty^z, so we need to add a
%minus sign!
pressure = -real(pressure);