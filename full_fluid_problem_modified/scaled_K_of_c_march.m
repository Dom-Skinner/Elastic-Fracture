%sets easy geometric parameters
n = 300;
t = round(n/2);
%xmax = 60;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;

%some values of lambda to try
%lambda = 0.056:0.0004:0.0592;
lambda = 0.0588:0.0004:0.0592;

coeff_data = zeros(4*n,length(lambda));

%error tolerance
tol = 10^(-7);

KI = zeros(1,length(lambda));

coeff_start = zeros(4*n,1);
coeff_start(n+1:4*n) = ones(3*n,1);

for i=1:length(lambda)
    if i == 2
        coeff_start = coeff_data(:,1);
    elseif i > 2
        coeff_start = (coeff_data(:,i-1)-coeff_data(:,i-2))/...
            (lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*coeff_data(:,i-2)-lambda(i-2)*...
            coeff_data(:,i-1))/(lambda(i-1)-lambda(i-2));
    end
    [KI(i),coeff_new,~] = ...
        scaled_fixed_lambda_M_iteration(n,xmax,lambda(i),tol,coeff_start);
    coeff_data(:,i) = coeff_new;
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
        coeff_start = (coeff_data(:,i-1)-coeff_data(:,i-2))/...
            (lambda(i-1)-lambda(i-2))*lambda(i) + (lambda(i-1)*...
            coeff_data(:,i-2)-lambda(i-2)*coeff_data(:,i-1))/...
            (lambda(i-1)-lambda(i-2));
    end    
    [KI(i),coeff_new,~] = ...
        scaled_fixed_lambda_M_iteration(n,xmax,lambda(i),tol,coeff_start);
    coeff_data(:,i) = hprime_new;
end