%sets easy geometric parameters
%n = 200;
t = round(n/2);
%xmax = 30;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;

%some values of lambda to try
%lambda = 0.056:0.0004:0.0592;
lambda = 0.056:0.0004:0.0588;

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

hprime_data = zeros(2*n,numel(lambda));
for k = 1:7
    hprime_data(1:n,k) = coeff_data(n+1:2*n,k) + x'.*coeff_data(1:n,k);
end

for k = 1:7
    hprime_data(n+1:2*n,k) = coeff_data(3*n+1:4*n,k) + x'.*coeff_data(2*n+1:3*n,k);
end