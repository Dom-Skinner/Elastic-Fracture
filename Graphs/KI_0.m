load n350x50-lambda

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