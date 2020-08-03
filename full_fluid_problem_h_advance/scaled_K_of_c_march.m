% For given parameters (lambda, L) and spacing x,z, iterate until static 
% solution is found. This is the basis for all of the other analysis.

% The values for x, z, determine somewhat the reliability/convergence
% properties of the whole solution. So they were created/experimented with
% in xend_test, and now we just load them.
clear
load n995x846-late
r = 180;
% These are the values of lambda and L used
%L = 0.8:0.2:3.6;
% 0.01;        early
% 0.2:0.2:0.6; mid
% 0.8:0.2:3.4; late
% 1:0.6:3.4 fixed
%lambda = 0:0.005:0.085;
% 0:0.008:0.058;    early
% 0.04:0.005:0.065; mid
% 0:0.005:0.085;    late
% 0.07:0.005:0.09;    fixedd


hprime_data = zeros(2*n+r,length(lambda),length(L));
tol = 10^(-8); %error tolerance
KI = zeros(length(lambda),length(L));
KII = zeros(length(lambda),length(L));

% Solve for each L, lambda combination
for j = 1:length(L)

[v,w] = spacing(x,z,L(j),r);

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
    [KI(i,j),KII(i,j),hprime_new,~] = ...
        scaled_fixed_lambda_M_iteration(x,z,v,w,t,x(end),lambda(i),tol,hprime_start);
    hprime_data(:,i,j) = hprime_new;
end


%{
s = 0.138673;
u = 4-6*s;
%%{
p1 = polyfit(KI(end-2:end,j).^u,lambda(end-2:end)',1);
l0(j) = p1(2);
p1 = polyfit(lambda(end-1:end)',KII(end-1:end,j),1);
p2 = polyfit(lambda(end-2:end)',KII(end-2:end,j),2);
KII0(j) = p2(1)*l0(j)*l0(j)+p2(2)*l0(j)+p2(3);
%}
%er = abs(p2(1)*l0(j)*l0(j)+p2(2)*l0(j)+p2(3) - (p1(1)*l0(j)+p1(2)));
% line above if you want some idea of the extrapolation error.

end

% plot the results for quick check:
subplot(1,3,1)
hold on
plot(L,KII.^2,'o-')
xlabel('\lambda_0')
ylabel('K_{II}')
subplot(1,3,2)
hold on
plot(lambda,KI,'o-')
xlabel('L')
ylabel('\lambda_0')
subplot(1,3,3)
hold on
plot(L,KII,'o-')
xlabel('L')
ylabel('K_{II}')