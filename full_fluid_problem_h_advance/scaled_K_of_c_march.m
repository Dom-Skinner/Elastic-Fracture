%sets easy geometric parameters
%%{
n = 350;
t = 170;
r = 50;
xmax = 50;

x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;


L = [0.6:0.3:2.4, 4:3:70];%[0.1,0.3,0.5,1, 1.5,2.1,2.6];



%some values of lambda to try
lambda =0.04:0.02:0.08;
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


%%{
s = 0.138673;
u = 4-6*s;
p1 = polyfit(KI(end-2:end,j).^u,lambda(end-2:end)',1);
l0(j) = p1(2);
p2 = polyfit(lambda(end-2:end)',KII(end-2:end,j),2);
KII0(j) = p2(1)*l0(j)*l0(j)+p2(2)*l0(j)+p2(3);

subplot(1,3,1)
plot(l0,KII0.^2,'o-')
xlabel('\lambda_0')
ylabel('K_{II}')
subplot(1,3,2)
plot(L,l0,'o-')
xlabel('L')
ylabel('\lambda_0')
subplot(1,3,3)
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
end