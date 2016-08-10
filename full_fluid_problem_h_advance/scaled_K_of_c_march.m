%sets easy geometric parameters
%%{
n = 250;
t = 120;
r = 20;
xmax = 50;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;


L = 2;
v = [tan((0:r-1)*atan(sqrt(L))/(r-1)).^2-L - 0.2 , x] + 2.2;
w = [tan((0.5:1:r-1.5)*atan(sqrt(L))/(r-1)).^2 - L - 0.2, 0.5*v(r) z]+2.2;





%some values of lambda to try
lambda =0.05:0.002:0.0588;
%}

%lambda = 0.0584:0.0004:0.0588;

hprime_data = zeros(2*n+r,length(lambda));

%error tolerance
tol = 10^(-8);

KI = zeros(1,length(lambda));

hprime_start = zeros(2*n+r,1);
hprime_start(1:n+r) = ones(n+r,1);
hprime_start(n+r+1:2*n+r) = x' + 1;
for i=1:length(lambda)
    if i == 2
        hprime_start = hprime_data(:,1);
    elseif i > 2
        hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/...
            (lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*hprime_data(:,i-2)-lambda(i-2)*...
            hprime_data(:,i-1))/(lambda(i-1)-lambda(i-2));
    end
    [KI(i),hprime_new,~] = ...
        scaled_fixed_lambda_M_iteration(x,z,v,w,t,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end

i = length(lambda); 