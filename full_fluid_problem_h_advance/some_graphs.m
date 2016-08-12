clf



%{
load n350x50-lambda

for j = 1:numel(lambda)
subplot(1,2,1)
plot(L, KI(j,:),'+-');
xlabel('L');
ylabel('K_I');
hold on
subplot(1,2,2)
plot(L, KII(j,:),'o-');
xlabel('L');
ylabel('K_{II}');
hold on
end
legend({num2str(lambda(1)),num2str(lambda(2)),num2str(lambda(3))})
%}
%subplot(1,3,1)

%{
load n350x50-lambda
for j = 1:numel(lambda)

subplot(1,2,1)
f1=KI(j,:).*(0.09982-0.0000651*(3.075-L).^6-lambda(j)).^(-1/u);
plot(L,f1,L,0.67-0.0001*(3.1-L).^6,'o-')
hold on
xlabel('L')
ylabel('K_{I}')
subplot(1,2,2)
plot(L,KII(j,:),'o-')
hold on
xlabel('L')
ylabel('K_{II}')
end
%}

%{
load n350x50-lambda

subplot(1,3,1)
plot(l0,KII0.^2,l0, 0.835*(0.09982-l0),'o-')
xlabel('\lambda_0')
ylabel('K_{II}')
subplot(1,3,2)
plot(L,l0,L,0.09982 - 0.0000651*(3.1-L).^6,'o-')
xlabel('L')
ylabel('\lambda_0')
subplot(1,3,3)
plot(L,KII0,L,0.0074*(3.1-L).^3,'o-')
xlabel('L')
ylabel('K_{II}')
%}
%{
for j = 3:numel(L)
plot(log(0.09982-lambda),log(KI(:,j)),log(0.09982-lambda),-0.39+0.32*log(0.09982-lambda),'o-')
hold on
end
%}