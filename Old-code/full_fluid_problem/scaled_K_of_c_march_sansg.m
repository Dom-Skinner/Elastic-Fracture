%sets easy geometric parameters
n = 100;
t = round(n/2);
xmax = 20;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;

%some values of lambda to try
lambda = 0:0.01:0.05;

hprime_data = zeros(n,length(lambda));

%error tolerance
tol = 10^(-8);

KI = zeros(1,length(lambda));

hprime_start = zeros(n,1);
hprime_start(1:n) = x' + 1;

for i=1:length(lambda)
    if i == 2
        hprime_start = hprime_data(:,1);
    elseif i > 2
        hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/(lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*hprime_data(:,i-2)-lambda(i-2)*hprime_data(:,i-1))/(lambda(i-1)-lambda(i-2));
    end

    [KI(i),hprime_new] = scaled_fixed_lambda_M_iteration_sansg(n,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end

i = length(lambda);

while 1
    i = i + 1;
    lambda(i) = (lambda(i-1)-lambda(i-2))/(KI(i-1)^3-KI(i-2)^3)*(KI(i-1)^3)/2 ...
        + (KI(i-1)^3*lambda(i-2)-KI(i-2)^3*lambda(i-1))/(KI(i-1)^3-KI(i-2)^3);
    hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/(lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*hprime_data(:,i-2)-lambda(i-2)*hprime_data(:,i-1))/(lambda(i-1)-lambda(i-2));
        
    [KI(i),hprime_new] = scaled_fixed_lambda_M_iteration_sansg(n,xmax,lambda(i),tol,hprime_start);
    
    hprime_data(:,i) = hprime_new;
end
