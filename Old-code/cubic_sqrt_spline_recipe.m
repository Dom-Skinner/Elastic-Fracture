function mixcoeffs = cubic_sqrt_spline_recipe(x, y, y1ddash, ynddash)



%splines with these values of x, y such that from x1 to xt, have
%representation 1/sqrt(x) (ax^3+bx^2+cx+d), and representation
%ax^3+bx^2+cx+d from xt+1 to xn.
%here, y(t) represents the value actually of sqrt(x)*actual function. this
%means that one of the splines must be adjusted accordingly.

n = length(x);
t = round(n/2);

%for splining on the square root side
x_sqrt = x(1:t);
y_sqrt = y(1:t);

%for splining on the pure side
x_pure = x(t:n);
y_pure = y(t:n);

x_join = x(t);

%modifies y_pure to match y_sqrt at xt
y_pure(1) = 1/(sqrt(x_join))*y(t);

%critically, we need to match up yt'' and yt'. This involves varying
%Ytddash_sqrt, Ytddash_pure, and matching them up them and Ytdash_sqrt,
%Ytdash_pure.

dd_sqrt_1 = 0;
dd_sqrt_2 = 1;
sqrt_coeffs_1 = cubic_spline_recipe(x_sqrt,y_sqrt,y1ddash,dd_sqrt_1);
sqrt_coeffs_2 = cubic_spline_recipe(x_sqrt,y_sqrt,y1ddash,dd_sqrt_2);

d_sqrt_1 = [(5/2)*x_join^(3/2), (3/2)*x_join^(1/2), (1/2)*x_join^(-1/2), ...
    (-1/2)*x_join^(-3/2)]*sqrt_coeffs_1(4*(t-1)-3:4*(t-1))';
d_sqrt_2 = [(5/2)*x_join^(3/2), (3/2)*x_join^(1/2), (1/2)*x_join^(-1/2), ...
    (-1/2)*x_join^(-3/2)]*sqrt_coeffs_2(4*(t-1)-3:4*(t-1))';

dd_sqrt_actual_1 = [(15/4)*x_join^(1/2), (3/4)*x_join^(-1/2), (-1/4)*x_join^(-3/2), ...
    (3/4)*x_join^(-5/2)]*sqrt_coeffs_1(4*(t-1)-3:4*(t-1))';
dd_sqrt_actual_2 = [(15/4)*x_join^(1/2), (3/4)*x_join^(-1/2), (-1/4)*x_join^(-3/2), ...
    (3/4)*x_join^(-5/2)]*sqrt_coeffs_2(4*(t-1)-3:4*(t-1))';

%gets actual coefficients for the linear variation, i.e.
A_dd_sqrt = (dd_sqrt_actual_2 - dd_sqrt_actual_1) ...
    /(dd_sqrt_2 - dd_sqrt_1);
B_dd_sqrt = (dd_sqrt_2*dd_sqrt_actual_1 - dd_sqrt_1*dd_sqrt_actual_2)...
    /(dd_sqrt_2 - dd_sqrt_1);
A_d_sqrt = (d_sqrt_2 - d_sqrt_1) ...
    /(dd_sqrt_2 - dd_sqrt_1);
B_d_sqrt = (dd_sqrt_2*d_sqrt_1 - dd_sqrt_1*d_sqrt_2)...
    /(dd_sqrt_2 - dd_sqrt_1);

%does the exact same but for the pure
dd_pure_1 = 0;
dd_pure_2 = 1;
pure_coeffs_1 = cubic_spline_recipe(x_pure,y_pure,dd_pure_1,ynddash);
pure_coeffs_2 = cubic_spline_recipe(x_pure,y_pure,dd_pure_2,ynddash);

d_pure_1 = [3*x_join^2, 2*x_join, 1, 0]*pure_coeffs_1(1:4)';
d_pure_2 = [3*x_join^2, 2*x_join, 1, 0]*pure_coeffs_2(1:4)';

A_dd_pure = 1;
B_dd_pure = 0;
A_d_pure = (d_pure_2 - d_pure_1) ...
    /(dd_pure_2 - dd_pure_1);
B_d_pure = (dd_pure_2*d_pure_1 - dd_pure_1*d_pure_2)...
    /(dd_pure_2 - dd_pure_1);

%solves the actual equations!!!
parameter_vector = [A_dd_sqrt, -A_dd_pure; A_d_sqrt, -A_d_pure]\[B_dd_pure - B_dd_sqrt; B_d_pure - B_d_sqrt];
Ytddash_sqrt = parameter_vector(1);
Ytddash_pure = parameter_vector(2);

mixcoeffs = zeros(1,4*(n-1),1);
mixcoeffs(1:4*(t-1)) = cubic_spline_recipe(x_sqrt,y_sqrt,y1ddash,Ytddash_sqrt);
mixcoeffs(4*(t-1)+1:4*(n-1)) = cubic_spline_recipe(x_pure,y_pure,Ytddash_pure,ynddash);
