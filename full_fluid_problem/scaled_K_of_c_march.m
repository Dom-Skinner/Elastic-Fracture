%sets easy geometric parameters
%n = 400;
t = round(n/2);
%xmax = 50;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;

%some values of lambda to try
lambda = 0.056:0.0004:0.0592;

hprime_data = zeros(2*n,length(lambda));

%error tolerance
tol = 10^(-8);

KI = zeros(1,length(lambda));

hprime_start = zeros(2*n,1);
hprime_start(1:n) = ones(n,1);
hprime_start(n+1:2*n) = x' + 1;
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
        scaled_fixed_lambda_M_iteration(n,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end

% Now we are done checking for the perscribed values of lambda, we will 
% now check as close as possible to lambda_0 until the program is manually
% terminated. Will be limited by numerical accuracy & the choice of n
i = length(lambda); 

while 0
    
    if KI(i) == -1 % If previous attempt failed, try again
        lambda(i) = [];
        return
        lambda(i) = lambda(i-1) + 0.5*(lambda(i)-lambda(i-1));
        
    else
        i=i+1;
        %
        lambda(i) = (lambda(i-1)-lambda(i-2))/(KI(i-1)^3-KI(i-2)^3)*...
            (KI(i-1)^3)/1.5 + (KI(i-1)^3*lambda(i-2)-KI(i-2)^3*...
            lambda(i-1))/(KI(i-1)^3-KI(i-2)^3);
        %
        hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/...
            (lambda(i-1)-lambda(i-2))*lambda(i) + (lambda(i-1)*...
            hprime_data(:,i-2)-lambda(i-2)*hprime_data(:,i-1))/...
            (lambda(i-1)-lambda(i-2));
    end    
    [KI(i),hprime_new,~] = ...
        scaled_fixed_lambda_M_iteration(n,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end