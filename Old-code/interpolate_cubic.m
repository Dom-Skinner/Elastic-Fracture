function matrix = interpolate_cubic(x)

%produces a matrix that interpolates the values of h', h'' to produce
%cubics.
%gets cubic representation up to panel t-1: i.e. uses
%x^(-1/2) * (ax^3 + bx^2 + cx + d)
%at panel t, has to ensure continuity of representation: cubic used here,
%but has to match up to later representation
%for rest of panels, uses representation (bx^2 + cx + d + e/x)
%doesn't do anything for crack tip: can extrapolate this from
%the first, and from the final panel.

%at infinity, uses another 9 points: for g = cx+d+e/x+f/x^2,
%and h = bx^2 + cx + d + e/x + f/x^2

%the input data is of the form: 1:n is g, n+1:2*n is g'. Use same one for
%h', obviously.
%however, for first t points of inpoint, input the g * x^(1/2)!!!
%this means the connecting panel is complicated.

%output of the form: in blocks of 4, for a, b, c, d respectively.

n = length(x);
t = round(n/2);

matrix = zeros(4*(n-1), 2*n);

%does the representation for first t-1 panels
for i = 1:t-1
    %need to solve a system of equations
    poly_matrix = [x(i)^3, x(i)^2, x(i)^1, x(i)^0; ... %for h(x(i))
        x(i+1)^3, x(i+1)^2, x(i+1)^1, x(i+1)^0; ... %for h(x(i+1)
        3*x(i)^2, 2*x(i), 1, 0; ... %for h'(x(i))
        3*x(i+1)^2, 2*x(i+1), 1, 0]; %for h'(x(i+1))
    matrix(4*(i-1)+1:4*i, [i, i+1, n+i, n+i+1]) = poly_matrix^(-1);
end

%does the connecting panel.
%use the x^(-1/2) * cubic rep here, but has to connect to other rep.
%solve system of equations
poly_matrix = [x(t)^3, x(t)^2, x(t)^1, 1; ... %for h(x(t))
    x(t+1)^(5/2), x(t+1)^(3/2), x(t)^(1/2), x(t)^(-1/2); ... %for h(x(t+1))
    3*x(t)^2, 2*x(t), 1, 0; ... %for h'(x(t))
    (5/2)*x(t+1)^(3/2), (3/2)*x(t+1)^(1/2), ...
    (1/2)*x(t+1)^(-1/2), (-1/2)*x(t+1)^(-3/2)]; %for h'(x(t+1))
matrix(4*(t-1)+1:4*t, [t, t+1, n+t, n+t+1]) = poly_matrix^(-1);

%does representation for second half of panels
for i = t+1:n-1
    %need to solve a system of equations
    poly_matrix = [x(i)^2, x(i), x(i)^(0), x(i)^(-1); ... %for h(x(i))
        x(i+1)^2, x(i+1), 1, x(i+1)^(-1); ... %for h(x(i+1)
        2*x(i), 1, 0, -x(i)^(-2); ... %for h'(x(i))
        2*x(i+1), 1, 0, -x(i+1)^(-2)]; %for h'(x(i+1))
    matrix(4*(i-1)+1:4*i, [i, i+1, n+i, n+i+1]) = poly_matrix^(-1);
end


