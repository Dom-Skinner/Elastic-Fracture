%we form a perturbation solution

%sets easy geometric parameters
n = 200;
t = round(n/2);
xmax = 40;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;

%a factor to scale h_0 by
h0prime = hprime_data(:,11);
lam0 = lambda(11);
%finds p0prime for this
h_values_matrix = hprime_to_h_values(x);
h0 = h_values_matrix*h0prime(n+1:2*n);
p0prime = zeros(n,1);
for j=1:n-1
        p0prime(j) = 4*lam0/(h0(j)+h0(j+1))^2;
end

scale = (lam0)^(-1/3);
%rescales
p0prime = scale*p0prime;
h0 = scale*h0;

%%%%
%make a 2nx2n matrix A for all the equations in h1'
%%%%
[kernel_matrix, interpolate_matrix] = pressure_shear_matrix(x,z);
shear_matrix = kernel_matrix(n:2*(n-1),:);
pressure_matrix = kernel_matrix(1:n-1,:);
A = zeros(2*n,2*n);

%need to exact p1' from p1
%we use the crude approximation p1'(x(n)) = 0
%p1prime_matrix = zeros((n-1),2*(n-1));
%for j=1:n-2
%        p1prime_matrix(j,j) = -1/(z(j+1)-z(j));
%        p1prime_matrix(j,j+1) = 1/(z(j+1)-z(j));
%end
%p1prime_matrix = p1prime_matrix*kernel_matrix;
%
%we need to produce a matrix for the eqn 2h0p0'h1 + h0^2p1' = 0, ie the
%fluid equation
%h_values_matrix_big = zeros(n,2*n);
%h_values_matrix_big(:,n+1:2*n) = h_values_matrix;
%fluid_matrix = zeros(n-1,2*n);
%for j=1:n-1
%    fluid_matrix(j,:) = 2*h0(j+1)*p0prime(j+1)*h_values_matrix_big(j+1,:) + ...
%        h0(j+1)^2*p1prime_matrix(j,:);
%end

%forms a matrix for p1, by integrating
%p1(z) = \int_z^{\infty} (2p0'/h0) h1
%takes in the coefficients of h1, and integrates
p1_matrix = zeros(n-1,3*n);
%runs over the panels. the final panel (the infinity panel) contributes
%nothing
for i=n-1:-1:1
    if i<n
        h0_midpanel = (h0(i)+h0(i+1))/2;
    else
        h0_midpanel = h0(i);
    end
    if i >= t
        for j=1:i
            if j < i
                p1_matrix(j,3*(i-1)+1) = (2*p0prime(i)/h0_midpanel)*(1/3)*(x(i+1)^3-x(i)^3);
                p1_matrix(j,3*(i-1)+2) = (2*p0prime(i)/h0_midpanel)*(1/2)*(x(i+1)^2-x(i)^2);
                p1_matrix(j,3*(i-1)+3) = (2*p0prime(i)/h0_midpanel)*(x(i+1)-x(i));
            else
                p1_matrix(j,3*(i-1)+1) = (2*p0prime(i)/h0_midpanel)*(1/3)*(x(i+1)^3-z(j)^3);
                p1_matrix(j,3*(i-1)+2) = (2*p0prime(i)/h0_midpanel)*(1/2)*(x(i+1)^2-z(j)^2);
                p1_matrix(j,3*(i-1)+3) = (2*p0prime(i)/h0_midpanel)*(x(i+1)-z(j));
            end
        end
    else
        for j=1:i
            if j < i
                p1_matrix(j,3*(i-1)+1) = (2*p0prime(i)/h0_midpanel)*(2/5)*(x(i+1)^(5/2)-x(i)^(5/2));
                p1_matrix(j,3*(i-1)+2) = (2*p0prime(i)/h0_midpanel)*(2/3)*(x(i+1)^(3/2)-x(i)^(3/2));
                p1_matrix(j,3*(i-1)+3) = (2*p0prime(i)/h0_midpanel)*(x(i+1)-x(i));
            else
                p1_matrix(j,3*(i-1)+1) = (2*p0prime(i)/h0_midpanel)*(2/5)*(x(i+1)^(5/2)-z(j)^(5/2));
                p1_matrix(j,3*(i-1)+2) = (2*p0prime(i)/h0_midpanel)*(2/3)*(x(i+1)^(3/2)-z(j)^(3/2));
                p1_matrix(j,3*(i-1)+3) = (2*p0prime(i)/h0_midpanel)*(x(i+1)-z(j));
            end      
        end
    end
end
%so far this matrix takes the coefficients of h1 as arguments. we convert
%it to take the coefficients of h1prime as arguments.
h_coefficient_matrix = hprime_to_h(x);
p1_matrix = p1_matrix*h_coefficient_matrix;
fluid_matrix = zeros(n-1,2*n);
fluid_matrix(1:n-1,n+1:2*n) = p1_matrix;


%makes the actual matrix A itself. first n-1 eqns are fluid.
A(1:n-1,:) = fluid_matrix - pressure_matrix;
A(n:2*(n-1),:) = shear_matrix;
%now next equations express that: h'' = 2g' at infinity
hdashdashrow = interpolate_matrix(3*n,:);
gdashrow = interpolate_matrix(2*n,:);
A(2*n-1,:) = hdashdashrow - 2*gdashrow;
%final equation expresses the K condition
A(2*n,:) = interpolate_matrix(3*n+1,:);


%the forcing equation
forcing = zeros(2*n,1);
forcing(2*n-1) = 0;
forcing(2*n) = 1;

h1prime = A\forcing;




