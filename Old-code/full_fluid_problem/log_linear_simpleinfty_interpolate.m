%a function to do interpolation for the coefficients of h'
%we ignore the final 5 values.

function interpolate_matrix = log_linear_simpleinfty_interpolate(x)

n = length(x);
t = round(n/2);

%creates a matrix for 4*n coefficients, incl 4 infinity coeffs
interpolate_matrix = zeros(4*n+5, 2*n);
for i = 1:t-1
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

%matching the two parts
%to find the gradients
interpolate_matrix(t,t) = -1/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(t,t+1) = 1/(x(t+1)-x(t));
interpolate_matrix(2*n+t,n+t) = -1/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(2*n+t,n+t+1) = 1/(x(t+1)-x(t));
%to find the constant terms
interpolate_matrix(n+t,t) = x(t+1)/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(n+t,t+1) = -x(t)/(x(t+1)-x(t));
interpolate_matrix(3*n+t,n+t) = x(t+1)/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(3*n+t,n+t+1) = -x(t)/(x(t+1)-x(t));

%for the purely linear part
for i = t+1:n-1
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

%linear term of g'
interpolate_matrix(n,:) = 0;
%constant term of g'
interpolate_matrix(2*n,n) = 0; %1;
%linear term of h'
interpolate_matrix(3*n,:) = 0; %interpolate_matrix(3*n-1,:);
%constant term of h'
interpolate_matrix(4*n, :) = 0; %interpolate_matrix(4*n-1,:);

interpolate_matrix(4*n+1:4*n+5,2*n+1:2*n+5) = eye(5);

return
end
