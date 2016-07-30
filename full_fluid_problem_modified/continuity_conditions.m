function continuity_mat = continuity_conditions(x)
% A function which returns a matrix representing the continuity conditions
% Should yield 0 = continuity_mat*coeffs
n = numel(x);
t = round(n/2);

continuity_mat = zeros(2*(n-1),4*n);

for k = 1:t-2
    continuity_mat(k,k)     =  x(k+1);
    continuity_mat(k,k+1)   = -x(k+1);
    continuity_mat(k,n+k)   =  1;
    continuity_mat(k,n+k+1) = -1;
    %
    % Same for h' coeffs
    continuity_mat(n-1+k,2*n+k)   =  x(k+1);
    continuity_mat(n-1+k,2*n+k+1) = -x(k+1);
    continuity_mat(n-1+k,3*n+k)   =  1;
    continuity_mat(n-1+k,3*n+k+1) = -1;
end

% Special continuity condition
continuity_mat(t-1,t-1)   = sqrt(x(t));
continuity_mat(t-1,t)     = -x(t);
continuity_mat(t-1,n+t-1) = 1/sqrt(x(t));
continuity_mat(t-1,n+t)   = -1;
% Almost identical for g'
continuity_mat(n-1+t-1,2*n+t-1) = sqrt(x(t));
continuity_mat(n-1+t-1,2*n+t)   = -x(t);
continuity_mat(n-1+t-1,3*n+t-1) = 1/sqrt(x(t));
continuity_mat(n-1+t-1,3*n+t)   = -1;

% and the rest

for k = t:n-1
    continuity_mat(k,k)     =  x(k+1);
    continuity_mat(k,k+1)   = -x(k+1);
    continuity_mat(k,n+k)   =  1;
    continuity_mat(k,n+k+1) = -1;
    %
    % Same for h' coeffs
    continuity_mat(n-1+k,2*n+k)   =  x(k+1);
    continuity_mat(n-1+k,2*n+k+1) = -x(k+1);
    continuity_mat(n-1+k,3*n+k)   =  1;
    continuity_mat(n-1+k,3*n+k+1) = -1;
end