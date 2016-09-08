%{
clear 
load n815x846
xmax = x(n);
%some values of lambda to try
lambda =0.0588* ( 1- (0:(1/12):1).^(2));
lambda(1:end) = lambda(end:-1:1);
scaled_K_of_c_march
%}
%return
load overall_fit

s = 0.138673;
l0 = 0.0591;
D = -0.00826;
K = 3*sqrt(2*pi)*KI;
K2 = 3*sqrt(2*pi)*hprime_data(1,:);
u = 4-6*s;
p3 = polyfit(K.^(1.5*u),lambda-l0-D*K.^u ,1);

er = max(abs(lambda-l0-D*K.^u-p3(1)*K.^(1.5*u)));
subplot(2,2,1)
plot(K,lambda,'*-',K,l0+D*K.^u+p3(1)*K.^(1.5*u),'o-')
fprintf('lambda = %.5e  %.5e K^u + %.5e K^(3*u/2)\n',l0,D,p3(1));
fprintf('With error %.5e \n\n', er);
xlabel('K_{I}')
ylabel('\lambda')
axis([0,2,0,0.06])

subplot(2,2,2)
p4 = polyfit(K2,lambda ,2);
plot(K2,lambda,K2,p4(1)*K2.^2 + p4(2)*K2 + p4(3),'o-')
fprintf('lambda =  %.5e KII^2 + %.5e KII + %.5e\n',p4(1),p4(2),p4(3));
fprintf('With error %.5e \n\n', max(abs(lambda-p4(1)*K2.^2-p4(2)*K2-p4(3))));
xlabel('K_{II}')
ylabel('\lambda')

fid = fopen('overall-fit-1.csv','w');
fprintf(fid,'KI,     KII,    lambda,  estimated,   estimated2\n');
for j = 1:numel(K)
    fprintf(fid, '%.5e,    %.5e,   %.5e,     %.5e,     %.5e\n', ...
        K(j), K2(j), lambda(j),l0+D*K(j).^u+p3(1)*K(j).^(1.5*u),...
        p4(1)*K2(j).^2 + p4(2)*K2(j) + p4(3));
end
fclose(fid);


clear

load n995x846-early
L1 = L; lam1 = l0; K2 = 3*sqrt(2*pi)*KII0;
load n995x846-mid
L1 = [L1 L]; lam1 = [lam1 l0]; K2 = [K2 3*sqrt(2*pi)*KII0];
load n995x846-late
L1 = [L1 L]; lam1 = [lam1 l0]; K2 = [K2 3*sqrt(2*pi)*KII0];

subplot(2,2,3)
p1 = polyfit(K2.^2,lam1 ,1);
fprintf('lambda = %.5e KII^2 + %.5e\n',p1(1),p1(2));
fprintf('With error %.5e \n\n', max(abs(lam1-p1(1)*K2.^2-p1(2))));
plot(K2,lam1,'o-',K2,p1(1)*K2.^2 +p1(2))
%axis([0,2,0.059,0.105])
ylabel('\lambda_0')
xlabel('K_{II}')




subplot(2,2,4)
hold on
p2 = polyfit(K2.^(1/5),L1 ,5);
plot(K2,L1,'o-',K2,real(p2(1)*K2 + p2(2)*K2.^(4/5) + p2(3)*K2.^(3/5) + ...
    p2(4)*K2.^(2/5)+p2(5)*K2.^(1/5) + p2(6)))
fprintf('L = %.5e K2 + %.5e *K2.^(4/5) + %.5e K2.^(3/5) %.5e  K2.^(2/5)+ %.5e K2.^(1/5) + %.5e\n' ...
    ,p2(1),p2(2),p2(3),p2(4),p2(5),p2(6))
fprintf('With error %.5e \n\n', max(abs(L1-real(p2(1)*K2 + ...
    p2(2)*K2.^(4/5) + p2(3)*K2.^(3/5) + ...
    p2(4)*K2.^(2/5)+p2(5)*K2.^(1/5) + p2(6)))));


xlabel('K_{II}')
ylabel('L')

fid = fopen('overall-fit-2.csv','w');
fprintf(fid,'KII,    lambda,  lambdaEst,  L,  LEst\n');
for j = 1:numel(K2)
    fprintf(fid, '%.5e,    %.5e,     %.5e,   %.5e,    %.5e\n', ...
        K2(j), lam1(j),p1(1)*K2(j).^2 +p1(2),...
        L1(j), p2(1)*K2(j) + p2(2)*K2(j).^(4/5) + p2(3)*K2(j).^(3/5) + ...
    p2(4)*K2(j).^(2/5)+p2(5)*K2(j).^(1/5) + p2(6));
end
fclose(fid);

