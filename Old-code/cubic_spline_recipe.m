%executes the cubic spline recipe found in Numerical Recipes.

%given points x1,x2,...,x(N),y1,y2,...,yN, y1'', yN'' finds the cubic
%spline through them
%outputs as coeff

function coeff = cubic_spline_recipe(x, y, y1ddash, yNddash)

N = length(x);

%forms and solves tridiagonal matrix.
mat_a = zeros(N,1);
mat_b = ones(N,1);
mat_c = zeros(N,1);
%fills in the diagonals
for j=2:N-1
    mat_a(j) = (x(j)-x(j-1))/6;
    mat_b(j) = (x(j+1)-x(j-1))/3;
    mat_c(j) = (x(j+1)-x(j))/6;
end
    
%finds the RHS
RHS_matrix = zeros(N,N+2);
RHS_matrix(1,N+1) = 1;
RHS_matrix(N,N+2) = 1;
for j = 2:N-1
    RHS_matrix(j,j-1) = 1/(x(j)-x(j-1));
    RHS_matrix(j,j) = -1/(x(j)-x(j-1)) -1/(x(j+1)-x(j));
    RHS_matrix(j,j+1) = 1/(x(j+1)-x(j));
end

yvector = zeros(N+2,1);
yvector(1:N) = y;
yvector(N+1) = y1ddash;
yvector(N+2) = yNddash;

RHS_vector = RHS_matrix*yvector;

%second derivative of y
yddash = tridiag(mat_a, mat_b, mat_c, RHS_vector, N);

coeff = zeros(1,4*(N-1));

%solves for the coefficients using y and yddash

for j = 1:N-1
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    
    %contribution of yj
    d = d + (x(j+1)/(x(j+1)-x(j)))*y(j);
    c = c - y(j)/(x(j+1)-x(j));
    
    %contribution of yj+1
    d = d - (x(j)/(x(j+1)-x(j)))*y(j+1);
    c = c + y(j+1)/(x(j+1)-x(j));
    
    %contribution of yj''
    d = d + yddash(j)*(6*(x(j+1)-x(j)))^(-1)*(x(j)*x(j+1)*(2*x(j+1)-x(j)));
    c = c - yddash(j)*(6*(x(j+1)-x(j)))^(-1)*(x(j)*x(j+1) + x(j+1)*(2*x(j+1)-x(j)) + x(j)*(2*x(j+1)-x(j)));
    b = b + yddash(j)*(6*(x(j+1)-x(j)))^(-1)*(x(j) + x(j+1) + 2*x(j+1) - x(j));
    a = a - yddash(j)*(6*(x(j+1)-x(j)))^(-1);
    
    %contribution of yj+1''
    d = d - yddash(j+1)*(6*(x(j+1)-x(j)))^(-1)*(x(j)*x(j+1)*(2*x(j)-x(j+1)));
    c = c + yddash(j+1)*(6*(x(j+1)-x(j)))^(-1)*(x(j)*x(j+1) + x(j+1)*(2*x(j)-x(j+1)) + x(j)*(2*x(j)-x(j+1)));
    b = b - yddash(j+1)*(6*(x(j+1)-x(j)))^(-1)*(x(j) + x(j+1) + 2*x(j) - x(j+1));
    a = a + yddash(j+1)*(6*(x(j+1)-x(j)))^(-1);
    
    coeff(4*(j-1)+4) = d;
    coeff(4*(j-1)+3) = c;
    coeff(4*(j-1)+2) = b;
    coeff(4*(j-1)+1) = a;
end
    
    
