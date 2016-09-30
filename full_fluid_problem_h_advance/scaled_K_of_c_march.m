%n = 350;
%t = 170;
r = 180;
xmax = x(n);



L = 0.8:0.2:3.4;
% 0.01;        early
% 0.2:0.2:0.6; mid
% 0.8:0.2:3.4; late
% 1:0.6:3.4 fixed



%some values of lambda to try
lambda = 0:0.005:0.085;
% 0:0.008:0.058;    early
% 0.04:0.005:0.065; mid
% 0:0.005:0.085;    late
% 0.07:0.005:0.09;    fixedd



hprime_data = zeros(2*n+r,length(lambda),length(L));

%error tolerance
tol = 10^(-8);

KI = zeros(length(lambda),length(L));
KII = zeros(length(lambda),length(L));
l0 = zeros(1,length(L));
KII0 = zeros(1,length(L));
hprime_l0 = zeros(t-1,length(L));

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
        scaled_fixed_lambda_M_iteration(x,z,v,w,t,xmax,lambda(i),tol,hprime_start);
    hprime_data(:,i,j) = hprime_new;
end


%%{
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
subplot(1,3,1)
hold on
plot(l0,KII0.^2,'o-')
xlabel('\lambda_0')
ylabel('K_{II}')
subplot(1,3,2)
hold on
plot(L,l0,'o-')
xlabel('L')
ylabel('\lambda_0')
subplot(1,3,3)
hold on
plot(L,KII0,'o-')
xlabel('L')
ylabel('K_{II}')

%{
for k = 1:t+r-1;
p3 = polyfit(KI(end-2:end,j).^2,hprime_data(k,end-2:end,j)',1);
hprime_l0(k,j) = p3(2);
end
hold on
plot(v(1:r+t-1)-L(j),v(1:r+t-1).^(-0.5)'.*hprime_l0(:,j));
%}
