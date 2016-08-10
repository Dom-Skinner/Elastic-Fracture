%takes in coefficients for h, stores as a 3*n vector, and gives out the
%values of p at the z values

function pressure = pprime_to_p(x, z, h_coeffs,t)

n = length(x);
infinity = 10^(10);

a = h_coeffs(1:3:3*n-2);
b = h_coeffs(2:3:3*n-1);
c = h_coeffs(3:3:3*n);


function ret = h_quad_int(x,a,b,c) 
ret = (-1).*(b.^2+(-4).*a.*c).^(-1).*((b+2.*a.*x).* ...
  (c+x.*(b+a.*x)).^(-1)+4.*a.*((-1).*b.^2+4.*a.*c).^(-1/2).* ...
  atan(((-1).*b.^2+4.*a.*c).^(-1/2).*(b+2.*a.*x)));
end
%does the h function

function ret = presum_sqrt_int(x,a,b,c,xi) 
    ret = ((4.*b.^2.* ...
  log(x.^(1/2)+(-1).*xi)+9.*a.*c.*log(x.^(1/2)+(-1).*xi).* ...
  xi).*(b+3.*a.*xi.^2).^(-1));
end

%roots_sqrt_int = @(a,b,c) cubic_root_custom(a,b,c);

function ret = h_sqrt_int(x,a,b,c) 
    root = cubic_root_custom(a,b,c);
    s1 = presum_sqrt_int(x,a,b,c,root(1,:));
    s2 = presum_sqrt_int(x,a,b,c,root(2,:));
    s3 = presum_sqrt_int(x,a,b,c,root(3,:));
    %
    ret = 2.*(4.*b.^3+27.*a.*c.^2).^(-1).*(c+x.^(1/2).* ...
    (b+a.*x)).^(-1).*(6.*b.*c+2.*b.^2.*x.^(1/2)+9.*a.*c.*x+(c+ ...
    x.^(1/2).*(b+a.*x)).*(s1+s2+s3));
end

%now we actually calculate the pressure at z(j)
h_sqrt_int_store_x  = h_sqrt_int(x(2:t),a(1:t-1)', b(1:t-1)',c(1:t-1)');
h_sqrt_int_store_z  = h_sqrt_int(z(1:t-1),a(1:t-1)', b(1:t-1)',c(1:t-1)');
h_sqrt_int_store_x0 = h_sqrt_int(x(1:t-1),a(1:t-1)', b(1:t-1)',c(1:t-1)');
%
%
h_quad_int_store_x  = h_quad_int([x(t+1:n) infinity],a(t:n)',b(t:n)',c(t:n)');
h_quad_int_store_z = h_quad_int(z(t:n-1),a(t:n-1)',b(t:n-1)',c(t:n-1)');
h_quad_int_store_x0  = h_quad_int(x(t:n) ,a(t:n)',b(t:n)',c(t:n)');

diag_x0_store = [h_sqrt_int_store_x - h_sqrt_int_store_x0, ...
    h_quad_int_store_x - h_quad_int_store_x0];
diag_z_store = [h_sqrt_int_store_x - h_sqrt_int_store_z, ...
    h_quad_int_store_x(1:end-1) - h_quad_int_store_z];

int_sum  = cumsum(diag_x0_store,'reverse');

pressure = diag_z_store +int_sum(2:end);


%oops we computed \int_z^\infty, not \int_infty^z, so we need to add a
%minus sign!
pressure = -real(pressure);
end