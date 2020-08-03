% This script looks for polynomial relationships between K,L and lambda

%sets easy geometric parameters
%{
n = 815;
t = 350;
r = 80;
xmax = 50;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;


L = 1:0.1:2;%[0.1,0.3,0.5,1, 1.5,2.1,2.6];



%some values of lambda to try
lambda =0.06:0.01:0.08;
%}

%lambda = 0.0584:0.0004:0.0588;

hprime_data = zeros(2*n+r,length(lambda),length(L));

%error tolerance
tol = 10^(-8);

KI = zeros(length(lambda),length(L));
KII = zeros(length(lambda),length(L));
l0 = zeros(1,length(L));
KII0 = zeros(1,length(L));
hprime_l0 = zeros(t-1,length(L));

for j = 1:length(L)

v = [tan((0:r)*atan(sqrt(L(j)))/r).^2 , x+L(j)] ;
w = [tan((0.5:1:r-0.5)*atan(sqrt(L(j)))/r).^2 - L(j) , 0.5*v(r) z]+L(j);
v(r+1)=[];
w(r+1)=[];

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


s = 0.138673;
u = 4-6*s;
p1 = polyfit(KI(end-2:end,j).^u,lambda(end-2:end)',1);
l0(j) = p1(2);
p2 = polyfit(lambda(end-2:end)',KII(end-2:end,j),2);
KII0(j) = p2(1)*l0(j)*l0(j)+p2(2)*l0(j)+p2(3);
end


p3 = polyfit(KII0.^2,l0,1);
l00 = p3(2);
c1 = p3(1);
% If in doubt about accuracy:
% plot(KII0.^2,l0,'o-',KII0.^2,KII0.^2 *c1+l00,'+-')
p4 = polyfit((l00-l0).^(1/6),L,1);
L0 = p4(2);
c2 = - (-1/p4(1))^6;
% If in doubt about accuracy:
% plot(L,l0,L,l00 +c2*(L0-L).^6)
c3 = sqrt(abs(c2/c1));

%plot(L,(lambda(1) - l00 - c2*(L0-L).^6)./(KI(1,:).^u))
c4 = mean((lambda(:) - l00 - c2*(L0-L(end)).^6)./(KI(:,end).^u));
% A bit of a rough approximation...