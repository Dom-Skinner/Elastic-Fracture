% For given parameter lambda and spacing x,z, iterate until static 
% solution is found. This is the basis for all of the other analysis.
% The result is the toughness constant KI, as well as the profile hprime,
% which can be used to find, e.g. the pressure using functions found in 
% other files.

clear

% The values for x, z, determine somewhat the reliability/convergence
% properties of the whole solution. We saved the best values, created by
% xend_test, but here we just create our own.
n = 1000;
t = 500;
xmax = 50;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
%some values of lambda to try
lambda = 0:0.012:0.0588;

hprime_data = zeros(2*n,length(lambda));

%error tolerance
tol = 10^(-8);

KI = zeros(1,length(lambda));

hprime_start = zeros(2*n,1);
hprime_start(1:n) = ones(n,1);
hprime_start(n+1:2*n) = x' + 1;
% Iterate through all lambda values, with updated hprime_start
for i=1:length(lambda)
    if i == 2
        hprime_start = hprime_data(:,1);
    elseif i > 2
        % To start the iteration, do an extrapolation of the previous
        % values. Not strictly necessary but will reduce the number of
        % iterations needed.
        hprime_start = (hprime_data(:,i-1)-hprime_data(:,i-2))/...
            (lambda(i-1)-lambda(i-2))*lambda(i) ...
            + (lambda(i-1)*hprime_data(:,i-2)-lambda(i-2)*...
            hprime_data(:,i-1))/(lambda(i-1)-lambda(i-2));
    end
    [KI(i),hprime_new,~] = ...
        scaled_fixed_lambda_M_iteration(x,z,n,t,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i) = hprime_new;
end



% plot the results for quick check:
plot(KI,lambda,'o-')
xlabel('\lambda_0')
ylabel('K_{I}')
