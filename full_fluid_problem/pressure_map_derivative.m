function deriv_pressure = pressure_map_derivative(x,z,h_coeffs,h_coefficient_matrix,t)
% WARNING: as a side effect of the optmized/vectorized nature of this code, 
% together with the fact that many of the mathematical formulas were found 
% symbolically, this code is extremely unreadable. It has been thoroughly
% tested though.

%returns the (n-1)xn matrix for the derivative of the map that sounds h' to
%pressure

%we first vary the 3n coefficients of h, to get a (n-1)x3n matrix. This is
%then multiplied by the map that takes h' to coefficients of h.
%h_coeffs = h_coefficient_matrix*hprime_old(n+1:2*n);
n = length(x);
infinity = 10^(10);

%the coefficients of h
a = h_coeffs(1:3:3*n-2);
b = h_coeffs(2:3:3*n-1);
c = h_coeffs(3:3:3*n);

%we first vary a. varying the first t-1 values of a variation require
%integrating -2/h^3 against x^(3/2) in panel t. The next require
%integrating -2/h^3 against x^2.

function ret = dh_a_sqrtrepn(x,a,b,c) 
ret = 2*(((2*b.^3 .*c.*sqrt(x) + b.^4.*x + b.^2.*c.*(c + a.*x) + ...
    ...
    8*a.*c.^2.*(c + 2*a.*x) + 2*a.*b.*c.*sqrt(x).*(5*c + 3*a.*x)) ... 
    ...
    ./(a.*(b.^2 - 4*a.*c).^2.*(c + b.*sqrt(x) + a.*x).^2)) ...
    ...
    + (12*b.*c.*atan((b + 2*a.*sqrt(x))./sqrt(-b.^2 + 4*a.*c))) ...
    ...
    ./(-b.^2 + 4*a.*c).^(5/2) );
end
%
%
    function ret = dh_a_quadrepn(x,a,b,c)
        ret = (-1).*(a.^(-1).*((-1).*b.^2+4.*a.*c).^(-1).*(b.*c+b.^2.*x+ ...
...
  (-2).*a.*c.*x).*(c+x.*(b+a.*x)).^(-2)+a.^(-1).*(b.^2+(-4).*a.*c) ...
...
  .^(-2).*(b.^2+2.*a.*c).*(b+2.*a.*x).*(c+x.*(b+a.*x)).^(-1)+4.*( ...
...
  b.^2+2.*a.*c).*((-1).*b.^2+4.*a.*c).^(-5/2).*atan(((-1).*b.^2+4.* ...
...
  a.*c).^(-1/2).*(b+2.*a.*x))); 
    end

function ret = dh_b_sqrtrepn(x,a,b,c) 

ret = -2*((b.*c + b.^2.*sqrt(x) - 2*a.*c.*sqrt(x))./(a.*(-b.^2 + 4*a.*c) ...
    ...
    .*(c + b.*sqrt(x) + a.*x).^2) + ((b.^2 + 2*a.*c) ...
    ...
    .*(b + 2*a.*sqrt(x)))./(a.*(b.^2 - 4*a.*c).^2.*(c + b.*sqrt(x) + ...
    ...
    a.*x)) + (4*(b.^2 + 2*a.*c).*atan((b + 2*a.*sqrt(x)) ...
    ...
    ./sqrt(-b.^2 + 4*a.*c)))./(-b.^2 + 4*a.*c).^(5/2));
end

    function ret = dh_b_quadrepn(x,a,b,c) 
        ret = (-1).*(b.^2+(-4).*a.*c).^(-2).*((b.^2+(-4).*a.*c).*(2.*c+ ...
...
  b.*x).*(c+x.*(b+a.*x)).^(-2)+(-3).*b.*(b+2.*a.*x).*(c+x.*(b+a.*x)) ...
...
  .^(-1)+(-12).*a.*b.*((-1).*b.^2+4.*a.*c).^(-1/2).*atan(((-1).* ...
 ...
 b.^2+4.*a.*c).^(-1/2).*(b+2.*a.*x)));
    end

function ret = dh_c_sqrtrepn(x,a,b,c) 
ret = 2*(((8*a.*c.^2 + 2*b.^3.*sqrt(x) + 2*a.*b.*sqrt(x) ...
    .*(5.*c + 3*a.*x) + b.^2.*(c + 9*a.*x))./((b.^2 - 4*a.*c).^2 ...
    .*(c + b.*sqrt(x) + a.*x).^2)) + (12*a.*b.*atan( ...
    (b + 2*a.*sqrt(x))./sqrt(-b.^2 + 4*a.*c)))./(-b.^2 + 4*a.*c).^(5/2));
end

    function ret = dh_c_quadrepn(x,a,b,c) 
        ret = (-1).*(b.^2+(-4).*a.*c).^(-2).*((-1).*(b+2.*a.*x).*(c+x.*( ...
...
  b+a.*x)).^(-2).*(b.^2+(-6).*a.*b.*x+(-2).*a.*(5.*c+3.*a.*x.^2))+ ...
...
  24.*a.^2.*((-1).*b.^2+4.*a.*c).^(-1/2).*atan(((-1).*b.^2+4.*a.*c) ...
...
  .^(-1/2).*(b+2.*a.*x)));
    end


%now produces a (n-1)x3n matrix with all these as entries

dp_h = zeros(n-1,3*n);

% All code now vectorised for speed with the downside that its almost 
% impossible to read
dh_a_sqrtrepn_store_x  = dh_a_sqrtrepn(x(2:t),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_a_sqrtrepn_store_x0 = dh_a_sqrtrepn(x(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_a_sqrtrepn_store_z  = dh_a_sqrtrepn(z(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
%
dh_b_sqrtrepn_store_x  = dh_b_sqrtrepn(x(2:t),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_b_sqrtrepn_store_x0 = dh_b_sqrtrepn(x(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_b_sqrtrepn_store_z  = dh_b_sqrtrepn(z(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
%
dh_c_sqrtrepn_store_x  = dh_c_sqrtrepn(x(2:t),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_c_sqrtrepn_store_x0 = dh_c_sqrtrepn(x(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
dh_c_sqrtrepn_store_z  = dh_c_sqrtrepn(z(1:t-1),a(1:t-1)',b(1:t-1)',c(1:t-1)');
%
%
dh_a_quadrepn_store_x  = dh_a_quadrepn([x(t+1:n) infinity],a(t:n)',b(t:n)',c(t:n)');
dh_a_quadrepn_store_x0 = dh_a_quadrepn(x(t:n) ,a(t:n)',b(t:n)',c(t:n)');
dh_a_quadrepn_store_z  = dh_a_quadrepn(z(t:n-1),a(t:n-1)',b(t:n-1)',c(t:n-1)');
%
dh_b_quadrepn_store_x  = dh_b_quadrepn([x(t+1:n) infinity],a(t:n)',b(t:n)',c(t:n)');
dh_b_quadrepn_store_x0 = dh_b_quadrepn(x(t:n) ,a(t:n)',b(t:n)',c(t:n)');
dh_b_quadrepn_store_z  = dh_b_quadrepn(z(t:n-1),a(t:n-1)',b(t:n-1)',c(t:n-1)');
%
dh_c_quadrepn_store_x  = dh_c_quadrepn([x(t+1:n) infinity],a(t:n)',b(t:n)',c(t:n)');
dh_c_quadrepn_store_x0 = dh_c_quadrepn(x(t:n),a(t:n)',b(t:n)',c(t:n)');
dh_c_quadrepn_store_z  = dh_c_quadrepn(z(t:n-1),a(t:n-1)',b(t:n-1)',c(t:n-1)');
%

diag_a_store = [dh_a_sqrtrepn_store_x - dh_a_sqrtrepn_store_z, ...
    dh_a_quadrepn_store_x(1:end-1) - dh_a_quadrepn_store_z];
diag_b_store = [dh_b_sqrtrepn_store_x - dh_b_sqrtrepn_store_z, ...
    dh_b_quadrepn_store_x(1:end-1) - dh_b_quadrepn_store_z];
diag_c_store = [dh_c_sqrtrepn_store_x - dh_c_sqrtrepn_store_z, ...
    dh_c_quadrepn_store_x(1:end-1) - dh_c_quadrepn_store_z];
%
diag_a_store0 = [dh_a_sqrtrepn_store_x - dh_a_sqrtrepn_store_x0, ...
        dh_a_quadrepn_store_x - dh_a_quadrepn_store_x0];
diag_b_store0 = [dh_b_sqrtrepn_store_x - dh_b_sqrtrepn_store_x0, ...
        dh_b_quadrepn_store_x - dh_b_quadrepn_store_x0];
diag_c_store0 = [dh_c_sqrtrepn_store_x - dh_c_sqrtrepn_store_x0, ...
        dh_c_quadrepn_store_x - dh_c_quadrepn_store_x0];


dp_h(:,1:3:end) = triu(repmat(diag_a_store0,n-1,1),1);
dp_h(:,2:3:end) = triu(repmat(diag_b_store0,n-1,1),1);
dp_h(:,3:3:end) = triu(repmat(diag_c_store0,n-1,1),1);

j = 1:n-1;
dp_h(3*(j-1)*(n-1)+j)     = diag_a_store;
dp_h((3*(j-1)+1)*(n-1)+j) = diag_b_store;
dp_h((3*(j-1)+2)*(n-1)+j) = diag_c_store;

%we computed it the wrong sign!
dp_h = -real(dp_h);

deriv_pressure = dp_h*h_coefficient_matrix;
end