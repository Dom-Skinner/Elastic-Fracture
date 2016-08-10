%a function to do interpolation for the coefficients of g'

function interpolate_matrix = linear_simpleinfty_interpolate_g(x,t,lambda)

n = length(x);

%creates a matrix for 4*n coefficients, incl 4 infinity coeffs
interpolate_matrix = zeros(2*n, n);
for i = 1:t-2
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
 
    %to find the constant terms
    interpolate_matrix(n+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n+i,i+1) = -x(i)/(x(i+1)-x(i));

    
end

%matching the two parts

%to find the gradients
interpolate_matrix(t-1,t-1) = -1/(x(t)-x(t-1));
interpolate_matrix(t-1,t) = sqrt(x(t))/(x(t)-x(t-1));


%to find the constant terms
interpolate_matrix(n+t-1,t-1) = x(t)/(x(t)-x(t-1));
interpolate_matrix(n+t-1,t) = -x(t-1)*sqrt(x(t))/(x(t)-x(t-1));

%for the purely linear part
for i = t:n-1
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpolate_matrix(n+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n+i,i+1) = -x(i)/(x(i+1)-x(i));
end

adjust_pen = 1-2*lambda/(3*x(n-1)) - 4*lambda^2 *log(x(n-1)) / x(n-1)^2;
adjust_end = 1-2*lambda/(3*x(n  )) - 4*lambda^2 *log(x(n  )) / x(n  )^2;
adjust_der = lambda/3/x(n)^2 + 4*lambda^2 *log(x(n))/x(n)^3;
%linear term of g'
interpolate_matrix(n,n) = 2*adjust_der/adjust_end ;
%constant term of g'
interpolate_matrix(2*n,n) = 1-x(n)*2*adjust_der/adjust_end ;

return
end
