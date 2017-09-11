%a function to do interpolation for the coefficients of h'

function interpolate_matrix = linear_simpleinfty_interpolate_h(x,t,lambda)

n = length(x);

%creates a matrix for 4*n coefficients, incl 4 infinity coeffs
interpolate_matrix = zeros(2*n, n);
for i = 1:t-2
    %to find the gradients

    interpolate_matrix(i,i) = -1/(sqrt(x(i+1))-sqrt(x(i)));
    interpolate_matrix(i,i+1) = 1/(sqrt(x(i+1))-sqrt(x(i)));
    %to find the constant terms
    interpolate_matrix(n+i,i) = sqrt(x(i+1))/(sqrt(x(i+1))-sqrt(x(i)));
    interpolate_matrix(n+i,i+1) = -sqrt(x(i))/(sqrt(x(i+1))-sqrt(x(i)));
    
end

%matching the two parts

%to find the gradients
interpolate_matrix(t-1,t-1) = -1/(sqrt(x(t))-sqrt(x(t-1)));
interpolate_matrix(t-1,t) = sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));

%to find the constant terms
interpolate_matrix(n+t-1,t-1) = sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));
interpolate_matrix(n+t-1,t) = -sqrt(x(t-1))*sqrt(x(t))/(sqrt(x(t))-sqrt(x(t-1)));

%for the purely linear part
for i = t:n-1
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpolate_matrix(n+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n+i,i+1) = -x(i)/(x(i+1)-x(i));
end

adjust_pen = 1-2*lambda/(pi*x(n-1));% - 4*lambda^2 *log(x(n-1)) / x(n-1)^2;
adjust_end = 1-2*lambda/(pi*x(n  ));% - 4*lambda^2 *log(x(n  )) / x(n  )^2;
%linear term of h'
interpolate_matrix(n,:) =  (adjust_end/adjust_pen)*interpolate_matrix(n-1,:);
%constant term of h'
interpolate_matrix(2*n, :) = interpolate_matrix(2*n-1,:) + ...
    x(n)*interpolate_matrix(n-1,:) -x(n)*interpolate_matrix(n,:);

return
end
