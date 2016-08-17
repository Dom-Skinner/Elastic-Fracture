load n350x50-lambda

for j = 1:3
subplot(2,2,1)
plot(L, KI(j,:),'+-');
xlabel('L');
ylabel('K_I');
hold on
subplot(2,2,2)
plot(L, KII(j,:),'o-');
xlabel('L');
ylabel('K_{II}');
hold on
end
legend({['\lambda=',num2str(lambda(1))],['\lambda=',num2str(lambda(2))],['\lambda=',num2str(lambda(3))]})

for j = 1:3
subplot(2,2,3)
plot(lambda, KI(:,j),'+-');
xlabel('\lambda');
ylabel('K_I');
hold on
subplot(2,2,4)
plot(lambda, KII(:,j),'o-');
xlabel('\lambda');
ylabel('K_{II}');
hold on
end
legend({['L=',num2str(L(1))],['L=',num2str(L(2))],['L=',num2str(L(3))]})