%takes in values y1, y2, ..., yn at points x1, x2, ..., xn along with a
%value of yn'' and y1'' = 0 and finds cubic coefficients in each panel

%again we use representation 1/sqrt(x) * (ax^3+bx^2+cx+d) up x(t+1)
%from x(t+1) to x(n), use ax^3+bx^2+cx+d

function [coeffs, forward_matrix] = cubic_spline_interpolation(x)

n = length(x);
t = round(n/2);
s = t+1; %the panel where the joining happens

%we should be able to, given y1,...,yn,yn'', y1''=0 determine
%a1,b1,c1,d1,...,an-1,bn-1,cn-1,dn-1 by inverting the right matrix.

%given such a1,....,dn-1, we evaluate at the n points, check 3*n-6
%continuinity conditions, and also y1'' and yn''

forward_matrix = zeros(4*(n-1),4*(n-1));
%rows of the forward matrix to give the y1,y2,...,yn. is adjusted for the
%two representations.

%rows for the value on the left: this one is simple
for i=1:n-1
    forward_matrix(i, 4*(i-1)+1:4*i) = [x(i)^3, x(i)^2, x(i), 1];
end

%for the value on the right: is harder due to mixed representation
for i=2:n
    if i == s
        forward_matrix(n-1+t, 4*(t-1)+1:4*t) = [x(s)^(5/2), x(s)^(3/2), x(s)^(1/2), x(s)^(-1/2)];
    else
        forward_matrix(n-1+(i-1), 4*(i-2)+1:4*(i-1)) = [x(i)^3, x(i)^2, x(i), 1];
    end
end

%check that y' is continuous: fiddly at s
for i=2:n-1
    if i == s
        forward_matrix(2*(n-1)+t, 4*t+1:4*s) = [3*x(s)^2, 2*x(s), 1, 0];
        forward_matrix(2*(n-1)+t, 4*(t-1)+1:4*t) = -[(5/2)*x(s)^(3/2), ...
            (3/2)*x(s)^(1/2), (1/2)*x(s)^(-1/2), (-1/2)*x(s)^(-3/2)];
    else
        forward_matrix(2*(n-1)+i-1, 4*(i-1)+1:4*i) = [3*x(i)^2, 2*x(i), 1, 0];
        forward_matrix(2*(n-1)+i-1, 4*(i-2)+1:4*(i-1)) = -[3*x(i)^2, 2*x(i), 1, 0];
    end
end

%check that y'' is continuous: fiddly at s
for i=2:n-1
    if i == s
        forward_matrix(2*(n-1)+n-2+t, 4*t+1:4*s) = [6*x(s), 2, 0, 0];
        forward_matrix(2*(n-1)+n-2+t, 4*(t-1)+1:4*t) = -[(15/4)*x(s)^(1/2), ...
            (3/4)*x(s)^(-1/2), (-1/4)*x(s)^(-3/2), (3/4)*x(s)^(-5/2)];
    else
        forward_matrix(2*(n-1)+n-2+i-1, 4*(i-1)+1:4*i) = [6*x(i), 2, 0, 0];
        forward_matrix(2*(n-1)+n-2+i-1, 4*(i-2)+1:4*(i-1)) = -[6*x(i), 2, 0, 0];
    end
end

%add condition on y''(x1)
forward_matrix(4*(n-1)-1,1:4) = [6*x(1), 2, 0, 0];

%add condition on y''(xn)
forward_matrix(4*(n-1),4*(n-1)-3:4*(n-1)) = [6*x(1), 2, 0, 0];

%prepares the data
preparation_matrix = zeros(4*(n-1),n+1);
preparation_matrix(1:n-1,1:n-1) = eye(n-1);
preparation_matrix(n:2*n-2, 2:n) = eye(n-1);
preparation_matrix(4*(n-1),n+1) = 1;

coeffs = forward_matrix^(-1)*preparation_matrix;


