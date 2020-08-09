% Look at catagory-edited.pdf to see what this file is for
clear
load K_lambda
K2 = 3*sqrt(2*pi)*hprime_data(1,:);
plot(K,K2,'o-')

hold on
plot([0,2.5],1.385*[1,1],'k')

plot([1.9322,2.5],1.5054*[1,1],'k')
plot(1.9322*[1,1],[1.5054,2.5],'k')

axis([0,2.5,1,2])
