%sets easy geometric parameters
n = 100;
t = round(n/2);
xmax = 20;

x = tan((0:n-1)*atan(xmax)/(n-1)).^2;

%some values of lambda to try
lambda(1:5)=0:5:20;

hprime_data = zeros(2*n,length(lambda));

%error tolerance
tol = 10^(-6);

KI = zeros(1,length(lambda));

hprime_start = zeros(2*n,1);
hprime_start(1:n) = ones(n,1);
hprime_start(n+1:2*n) = x' + 1;

for i=1:length(lambda)
    if i == 2
        hprime_start = hprime_data(:,1);
    elseif i > 2
        hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/(lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*hprime_data(:,i-2)-lambda(i-2)*hprime_data(:,i-1))/(lambda(i-1)-lambda(i-2));
    end
    [KI(i),hprime_new] = fixed_lambda_M_iteration(n,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end