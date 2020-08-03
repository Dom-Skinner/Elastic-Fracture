%a function to do interpolation for the coefficients of h'

function interpolate_matrix = linear_simpleinfty_interpolate(x,t,lambda)

n = length(x);

%creates a matrix for 4*n coefficients, incl 4 infinity coeffs
interpolate_matrix = zeros(4*n, 2*n);
for i = 1:t-2
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    interpolate_matrix(2*n+i,n+i) = -1/(sqrt(x(i+1))-sqrt(x(i)));
    interpolate_matrix(2*n+i,n+i+1) = 1/(sqrt(x(i+1))-sqrt(x(i)));
    %to find the constant terms
    interpolate_matrix(n+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n+i,i+1) = -x(i)/(x(i+1)-x(i));
    interpolate_matrix(3*n+i,n+i) = sqrt(x(i+1))/(sqrt(x(i+1))-sqrt(x(i)));
    interpolate_matrix(3*n+i,n+i+1) = -sqrt(x(i))/(sqrt(x(i+1))-sqrt(x(i)));
    
end

%matching the two parts

%to find the gradients
interpolate_matrix(t-1,t-1) = -1/(x(t)-x(t-1));
interpolate_matrix(t-1,t) = sqrt(x(t))/(x(t)-x(t-1));
interpolate_matrix(2*n+t-1,n+t-1) = -1/(sqrt(x(t))-sqrt(x(t-1)));
interpolate_matrix(2*n+t-1,n+t) = sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));

%to find the constant terms
interpolate_matrix(n+t-1,t-1) = x(t)/(x(t)-x(t-1));
interpolate_matrix(n+t-1,t) = -x(t-1)*sqrt(x(t))/(x(t)-x(t-1));
interpolate_matrix(3*n+t-1,n+t-1) = sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));
interpolate_matrix(3*n+t-1,n+t) = -sqrt(x(t-1))*sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));

%for the purely linear part
for i = t:n-1
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    interpolate_matrix(2*n+i,n+i) = -1/(x(i+1)-x(i));
    interpolate_matrix(2*n+i,n+i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpolate_matrix(n+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n+i,i+1) = -x(i)/(x(i+1)-x(i));
    interpolate_matrix(3*n+i,n+i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(3*n+i,n+i+1) = -x(i)/(x(i+1)-x(i));
end

adjust_pen = 1-2*lambda/(pi*x(n-1));
adjust_end = 1-2*lambda/(pi*x(n  ));
adjust_der = lambda/pi/x(n)^2;
%linear term of g'
interpolate_matrix(n,n) = 2*adjust_der/adjust_end ;
%constant term of g'
interpolate_matrix(2*n,n) = 1-x(n)*2*adjust_der/adjust_end ;
%linear term of h'
interpolate_matrix(3*n,:) =  (adjust_end/adjust_pen)*interpolate_matrix(3*n-1,:);
%constant term of h'
interpolate_matrix(4*n, :) = interpolate_matrix(4*n-1,:) + ...
    x(n)*interpolate_matrix(3*n-1,:) -x(n)*interpolate_matrix(3*n,:);

return
end
